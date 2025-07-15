import os, sys
from pathlib import Path
from statistics import mean
import glob
import pandas as pd
import numpy as np
import requests
from typing import Dict, Any, Optional


def calculate_delta_psi(sum_bg_file, bg_file, target_file, outdir):
    """
    Calculate delta PSI values and p-values for GTEx TCGA exon comparisons.
    
    Parameters:
    - sum_bg_file: path to the GTEx summarized background file (one psi column + grouping cols)
    - bg_file: path to the raw GTEx spliced exon file (multiple psi sample columns)
    - target_file: path to the TCGA spliced exon file (multiple psi sample columns)
    - outdir: directory where outputs will be saved
    
    Outputs:
    - {outdir}/{sample}-psi.txt: TSV with columns
      [ensembl_gene_id, symbol, chr, strand, exon_start, exon_end,
       psi_gtex, psi_tcga, delta_psi, pval]
    """
    # Helper: empirical CDF generator
    def get_ecdf(vals):
        v = np.asarray(vals, dtype=float)
        valid = ~np.isnan(v)
        n = valid.sum()
        min_valid = max(10, int(n * 0.1))
        if n < min_valid:
            return lambda x: np.nan
        diffs = np.abs(v[valid][:, None] - v[valid][None, :])
        i, j = np.tril_indices(n, k=-1)
        subtracts = diffs[i, j]
        if np.all(subtracts == 0):
            return lambda x: np.nan
        sorted_subs = np.sort(subtracts)
        m = len(sorted_subs)
        def ecdf(x):
            # fraction of background diffs <= x
            return np.searchsorted(sorted_subs, x, side='right') / m
        return ecdf

    # Ensure output directory exists
    os.makedirs(outdir, exist_ok=True)
    
    # Define grouping columns
    group_cols = ['ensembl_gene_id', 'symbol', 'chr', 'strand', 'exon_start', 'exon_end']

    # 1) Load and preprocess GTEx summary (one PSI column)
    sum_bg = pd.read_csv(sum_bg_file, sep='\t')
    sum_bg = sum_bg.groupby(group_cols, as_index=False).mean()
    psi_cols = [c for c in sum_bg.columns if c not in group_cols]
    if len(psi_cols) != 1:
        raise ValueError("Expected exactly one PSI column in summarized background file.")
    sum_bg = sum_bg.rename(columns={psi_cols[0]: 'psi_gtex'})

    # 2) Load and preprocess raw GTEx for ECDF
    bg_all = pd.read_csv(bg_file, sep='\t')
    bg_all = bg_all.drop(columns=['upstreamEE','downstreamES'], errors='ignore')
    bg_all = bg_all.groupby(['GeneID','geneSymbol','chr','strand','exonStart_0base','exonEnd'], as_index=False).mean()
    bg_all['chr'] = bg_all['chr'].astype(str).str.replace('chr','')
    bg_all['strand'] = bg_all['strand'].astype(str).map(lambda s: -1 if s=='-' else 1)
    bg_all = bg_all.rename(columns={
        'GeneID':'ensembl_gene_id',
        'geneSymbol':'symbol',
        'exonStart_0base':'exon_start',
        'exonEnd':'exon_end'
    })
    sample_cols_bg = [c for c in bg_all.columns if c not in group_cols]
    bg_all = bg_all[group_cols + sample_cols_bg]
    
    # Precompute ECDF functions per exon
    ecdf_funcs = [get_ecdf(row[sample_cols_bg].values) for _, row in bg_all.iterrows()]
    ecdf_df = bg_all[group_cols].copy()
    ecdf_df['ecdf_func'] = ecdf_funcs

    # 3) Load and preprocess TCGA data
    tcga_all = pd.read_csv(target_file, sep='\t')
    tcga_all = tcga_all.drop(columns=['upstreamEE','downstreamES'], errors='ignore')
    tcga_all = tcga_all.groupby(['AC','GeneName','chr','strand','exonStart','exonEnd'], as_index=False).mean()
    tcga_all['chr'] = tcga_all['chr'].astype(str).str.replace('chr','')
    tcga_all['strand'] = tcga_all['strand'].astype(str).map(lambda s: -1 if s=='-' else 1)
    tcga_all = tcga_all.rename(columns={
        'AC':'ensembl_gene_id',
        'GeneName':'symbol',
        'exonStart':'exon_start',
        'exonEnd':'exon_end'
    })
    sample_cols_tcga = [c for c in tcga_all.columns if c not in group_cols]
    tcga_all = tcga_all[group_cols + sample_cols_tcga]

    # 4) For each TCGA sample: compute delta PSI and p-values
    for sam in sample_cols_tcga:
        print(f"Processing sample {sam}...")
        tmp = tcga_all[group_cols + [sam]].copy()
        tmp = tmp.rename(columns={sam: 'psi_tcga'})
        
        # merge with GTEx summary
        joint = pd.merge(sum_bg, tmp, on=group_cols, how='inner')
        joint = joint.dropna(subset=['psi_gtex','psi_tcga'])
        joint['delta_psi'] = joint['psi_tcga'] - joint['psi_gtex']
        
        # merge with ECDF functions
        joint = pd.merge(joint, ecdf_df, on=group_cols, how='left')
        joint['abs_delta'] = np.abs(joint['delta_psi'])
        
        # evaluate CDF and compute p-value
        joint['cdf_val'] = joint.apply(
            lambda row: row['ecdf_func'](row['abs_delta']) if callable(row['ecdf_func']) else np.nan,
            axis=1
        )
        joint['pval'] = 1 - joint['cdf_val']
        
        # keep only numeric columns; drop the function pointers
        out_df = joint.drop(columns=['ecdf_func','abs_delta','cdf_val'])
        
        # write to TSV
        out_path = os.path.join(outdir, f"{sam}-psi.txt")
        out_df.to_csv(out_path, sep='\t', index=False)


def combine_spliced_exon(in_dir):

  def get_nums_list(s):
      '''
      given a list of string numbers
      with possible NaNs return
      all non-NaN values
      '''
      l = []

      for elm in s:
          if elm != 'NaN' and elm.strip() != "":
              l.append(float(elm))

      return l


  for fn in os.listdir(in_dir):
      if fn.endswith('combined_mean.txt'): # skip generated output files 
          continue
      print("processing file: " + fn)
      out_mean = open(in_dir + "/" + fn.replace('.txt', '') + '_combined_mean.txt', 'w')
      out_mean.write('ensembl_gene_id\tsymbol\tchr\tstrand\texon_start\texon_end\tpsi\n')

      header = True
      with open(in_dir + "/" + fn, 'r') as f:
          for line in f:
              if header:
                  header = False
                  continue
              words = line.strip().split('\t')
              strand = "1"
              if words[3] == '-':
                  strand = "-1"
              s = words[0] + "\t" + words[1] + "\t" + words[2].replace("chr", "") + "\t" + strand + "\t" + words[4] + '\t' + words[5] + "\t"
              nums = get_nums_list(words[8:])
              a = "NaN"
              if len(nums) > 1:
                  a = mean(nums)
              out_mean.write(s + str(a) + "\n")
      out_mean.close()



def preprocess_pooled(
    compare_file: str,
    out_psi_dir: str,
    *,
    tissue: Optional[str] = None,
    tissue_download_root: Optional[str] = None,
    normal_path: Optional[str] = None
) -> Dict[str, Any]:
    """
    Either download a GTEx splicing matrix for `tissue` or use your own `normal_path`,
    then combine spliced exons and compute delta PSI.

    Args:
        compare_file:   Path to your compare.txt output file.
        out_psi_dir:        Output directory containing PSI files that can be inputed to splitpea.

        tissue:         (optional) Tissue name from IRIS that will be dowloaded. 
                        Y. Pan, J.W. Phillips, B.D. Zhang, M. Noguchi, E. Kutschera, J. McLaughlin, P.A. Nesterenko, Z. Mao, N.J. Bangayan, R. Wang, W. Tran, H.T. Yang, Y. Wang, Y. Xu, M.B. Obusan, D. Cheng, A.H. Lee, K.E. Kadash-Edmondson, A. Champhekar, C. Puig-Saus, A. Ribas, R.M. Prins, C.S. Seet, G.M. Crooks, O.N. Witte, & Y. Xing, IRIS: Discovery of cancer immunotherapy targets arising from pre-mRNA alternative splicing, Proc. Natl. Acad. Sci. U.S.A. 120 (21) e2221116120, https://doi.org/10.1073/pnas.2221116120 (2023).

        tissue_download_root:  (required if `tissue` is set) Root directory under which to
                        create GTEx_<Tissue>/splicing_matrix/.

        normal_path:    (optional) Path to a pre-downloaded normal splicing matrix file.
    
    e.g:
        preprocess_pooled(
            tissue="uterus",
            download_root="/uterus",
            compare_file="splicing_matrix.SE.cov10.TCGA_UCS_T.txt",
            psi_dir="/UCS_uterus"
        )
    """
    if bool(tissue) == bool(normal_path):
        raise ValueError("You must provide exactly one of `tissue` or `normal_path`.")
    
    if normal_path:
        if not os.path.isfile(normal_path):
            raise FileNotFoundError(f"normal_path not found: {normal_path}")
        normal_file = normal_path
        work_dir    = os.path.dirname(normal_path)

    else:
        _ALLOWED_TISSUES = [
            "AdiposeTissue", "AdrenalGland", "Bladder", "Blood", "BloodVessel",
            "Brain", "Breast", "CervixUteri", "Colon", "Esophagus", "FallopianTube",
            "Heart", "Kidney", "Liver", "Lung", "Muscle", "Nerve", "Ovary",
            "Pancreas", "Pituitary", "Prostate", "SalivaryGland", "Skin",
            "SmallIntestine", "Spleen", "Stomach", "Testis", "Thyroid",
            "Uterus", "Vagina"
        ]
        slug_map = {
            t.lower().replace(" ", "").replace("_", ""): t
            for t in _ALLOWED_TISSUES
        }

        slug = tissue.lower().replace(" ", "").replace("_", "")
        if slug not in slug_map:
            raise ValueError(
                f"Unknown tissue '{tissue}'. Valid choices are: {', '.join(_ALLOWED_TISSUES)}"
            )
        suffix = slug_map[slug]

        if not tissue_download_root:
            raise ValueError("`tissue_download_root` must be set when using `tissue`.")

        remote_base = (
            "https://xinglabtrackhub.research.chop.edu"
            "/iris/IRIS_data.v2.0.0/db"
        )
        fname       = f"splicing_matrix.SE.cov10.GTEx_{suffix}.txt"
        url         = f"{remote_base}/GTEx_{suffix}/splicing_matrix/{fname}"

        work_dir    = os.path.join(tissue_download_root, f"GTEx_{suffix}", "splicing_matrix")
        os.makedirs(work_dir, exist_ok=True)
        normal_file = os.path.join(work_dir, fname)

        if not os.path.exists(normal_file):
            print(f"Downloading {url} to {normal_file}")
            print("Y. Pan, J.W. Phillips, B.D. Zhang, M. Noguchi, E. Kutschera, J. McLaughlin, P.A. Nesterenko, Z. Mao, N.J. Bangayan, R. Wang, W. Tran, H.T. Yang, Y. Wang, Y. Xu, M.B. Obusan, D. Cheng, A.H. Lee, K.E. Kadash-Edmondson, A. Champhekar, C. Puig-Saus, A. Ribas, R.M. Prins, C.S. Seet, G.M. Crooks, O.N. Witte, & Y. Xing, IRIS: Discovery of cancer immunotherapy targets arising from pre-mRNA alternative splicing, Proc. Natl. Acad. Sci. U.S.A. 120 (21) e2221116120, https://doi.org/10.1073/pnas.2221116120 (2023).")
            resp = requests.get(url, stream=True)
            resp.raise_for_status()
            with open(normal_file, "wb") as wf:
                for chunk in resp.iter_content(chunk_size=8192):
                    wf.write(chunk)

    combine_spliced_exon(work_dir)

    combined_list = glob.glob(os.path.join(work_dir, "*_combined_mean.txt"))
    if not combined_list:
        raise FileNotFoundError(f"No *_combined_mean.txt found in {work_dir}")
    combined_file = combined_list[0]

    print("Starting delta PSI calculation...")
    delta = calculate_delta_psi(
        combined_file,  
        normal_file,    
        compare_file,
        out_psi_dir
    )

    return out_psi_dir