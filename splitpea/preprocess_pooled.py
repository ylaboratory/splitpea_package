import os, sys
from pathlib import Path
from statistics import mean
import re
import glob
import pandas as pd
import numpy as np
import requests
from typing import Dict, Any, Optional

import importlib_resources as pkg_resources
from .src import reference
from .src import mouse_ref

def calculate_delta_psi(sum_bg_file, bg_file, target_file, outdir):
    """
    Calculate delta PSI values and p-values for background_samples compare_samples exon comparisons.
    
    Parameters:
    - sum_bg_file: path to the background_samples summarized background file (one psi column + grouping cols)
    - bg_file: path to the raw background_samples spliced exon file (multiple psi sample columns)
    - target_file: path to the compare_samples spliced exon file (multiple psi sample columns)
    - outdir: directory where outputs will be saved
    
    Outputs:
    - {outdir}/{sample}-psi.txt: TSV with columns
      [ensembl_gene_id, symbol, chr, strand, exon_start, exon_end,
       psi_background_samples, psi_compare_samples, delta_psi, pval]
    """
    os.makedirs(outdir, exist_ok=True)

    group_cols = ['ensembl_gene_id', 'symbol', 'chr', 'strand', 'exon_start', 'exon_end']

    def normalize_bg(df):
        df = df.rename(columns={'exonStart': 'exonStart_0base'})
        df = df.drop(columns=['upstreamEE', 'downstreamES'], errors='ignore')
        df = df.groupby(['GeneID','geneSymbol','chr','strand','exonStart_0base','exonEnd'], as_index=False).mean()
        df['chr'] = df['chr'].astype(str).str.replace('chr', '', regex=False)
        df['strand'] = df['strand'].astype(str).map(lambda s: -1 if s == '-' else 1)
        df = df.rename(columns={
            'GeneID': 'ensembl_gene_id',
            'geneSymbol': 'symbol',
            'exonStart_0base': 'exon_start',
            'exonEnd': 'exon_end'
        })
        return df

    def normalize_cmp(df):
        df = df.rename(columns={'GeneID': 'AC', 'geneSymbol': 'GeneName'})
        df = df.drop(columns=['upstreamEE', 'downstreamES'], errors='ignore')
        df = df.groupby(['AC','GeneName','chr','strand','exonStart','exonEnd'], as_index=False).mean()
        df['chr'] = df['chr'].astype(str).str.replace('chr', '', regex=False)
        df['strand'] = df['strand'].astype(str).map(lambda s: -1 if s == '-' else 1)
        df = df.rename(columns={
            'AC': 'ensembl_gene_id',
            'GeneName': 'symbol',
            'exonStart': 'exon_start',
            'exonEnd': 'exon_end'
        })
        return df

    def count_pairs_leq_diff(v_sorted, x):
        n = len(v_sorted)
        if n < 2:
            return 0
        j = 0
        total = 0
        for i in range(n):
            while j < n and v_sorted[j] - v_sorted[i] <= x:
                j += 1
            total += max(0, (j - i - 1))
        return total

    sum_bg = pd.read_csv(sum_bg_file, sep='\t')
    sum_bg = sum_bg.groupby(group_cols, as_index=False).mean()
    psi_cols = [c for c in sum_bg.columns if c not in group_cols]
    if len(psi_cols) != 1:
        raise ValueError("Expected exactly one PSI column in summarized background file.")
    sum_bg = sum_bg.rename(columns={psi_cols[0]: 'psi_background_samples'})

    bg_all = pd.read_csv(bg_file, sep='\t')
    bg_all = normalize_bg(bg_all)
    sample_cols_bg = [c for c in bg_all.columns if c not in group_cols]
    bg_all = bg_all[group_cols + sample_cols_bg]

    bg_num = bg_all[sample_cols_bg].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)

    ens_   = bg_all["ensembl_gene_id"].to_numpy()
    sym_   = bg_all["symbol"].to_numpy()
    chr_  = bg_all["chr"].to_numpy()
    strand_ = bg_all["strand"].astype("int64").to_numpy()
    exon_start_ = bg_all["exon_start"].astype("int64").to_numpy()
    exon_end_   = bg_all["exon_end"].astype("int64").to_numpy()

    keys = list(zip(ens_, sym_, chr_, strand_, exon_start_, exon_end_))

    bg_groups = {}
    for key, row_vals in zip(keys, bg_num):
        vals = row_vals[~np.isnan(row_vals)]
        n = vals.size
        min_valid = max(10, int(n * 0.1)) 
        if n >= min_valid:
            bg_groups[key] = np.sort(vals)

    compare = pd.read_csv(target_file, sep='\t')
    compare = normalize_cmp(compare)
    sample_cols_cmp = [c for c in compare.columns if c not in group_cols]
    compare = compare[group_cols + sample_cols_cmp]

    cols_out = group_cols + ['psi_background_samples', 'psi_compare_samples', 'delta_psi', 'pval']

    for sam in sample_cols_cmp:
        print(f"Processing sample {sam}...")
        tmp = compare[group_cols + [sam]].copy()
        tmp = tmp.rename(columns={sam: 'psi_compare_samples'})
        joint = pd.merge(sum_bg, tmp, on=group_cols, how='inner', sort=False)
        joint = joint.dropna(subset=['psi_background_samples', 'psi_compare_samples'])
        joint = joint.reset_index(drop=True)
        joint['delta_psi'] = joint['psi_compare_samples'] - joint['psi_background_samples']
        joint['abs_delta'] = np.abs(joint['delta_psi'])
        joint['_exon_key'] = joint[group_cols].apply(tuple, axis=1)
        pval_cache = {}
        pvals = np.full(len(joint), np.nan, dtype=float)
        for exon_key, idx_series in joint.groupby('_exon_key').groups.items():
            idx = np.fromiter(idx_series, dtype=int)
            v_sorted = bg_groups.get(exon_key, None)
            if v_sorted is None or len(v_sorted) < 2 or v_sorted[-1] == v_sorted[0]:
                continue
            n = len(v_sorted)
            m_pairs = n * (n - 1) // 2
            if m_pairs == 0:
                continue
            xvals = np.asarray(sorted(joint['abs_delta'].iloc[idx].dropna().unique()))
            if xvals.size == 0:
                continue
            cdf_map = {}
            for x in xvals:
                key2 = (exon_key, float(x))
                if key2 in pval_cache:
                    cdf = pval_cache[key2]
                else:
                    cnt = count_pairs_leq_diff(v_sorted, x)
                    cdf = cnt / m_pairs
                    pval_cache[key2] = cdf
                cdf_map[x] = cdf
            abs_deltas = joint['abs_delta'].iloc[idx].to_numpy()
            pv = np.array([1.0 - cdf_map.get(x, np.nan) for x in abs_deltas], dtype=float)
            pvals[idx] = pv
        joint['pval'] = pvals
        out_df = joint.drop(columns=['abs_delta', '_exon_key'])
        out_path = os.path.join(outdir, f"{sam}-psi.txt")
        out_df[cols_out].to_csv(out_path, sep='\t', index=False)



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


def read_sample_specific(
    in_file: str,
    map_path: str = None,
    species: str = "human",
    single_rMATS: bool = False,
    inclevel: int = 1
) -> pd.DataFrame:
    """
    Convert rMATS output into a unified PSI table.

    Parameters:
    - in_file: path to a directory of rMATS files, or (if single_rMATS=True) a single rMATS .txt file
    - map_path: optional Ensembl→symbol mapping TSV with columns [symbol, ensembl]
    - species: "human" or "mouse" (used if map_path not given)
    - single_rMATS: if True, treat `in_file` as a single rMATS file instead of a directory
    - inclevel: which rMATS inclusion-level field to use when parsing a single file; 1 or 2. Defaults to 1 if not specified.

    Outputs:
    - DataFrame columns:
      [GeneID, geneSymbol, chr, strand, exonStart, exonEnd,
       upstreamEE, downstreamES, sample1, sample2, ...]
    """

    def norm_id(x: str) -> str:
        if pd.isna(x):
            return x
        return str(x)

    def fix_chr(x: str):
        if pd.isna(x):
            return x
        s = str(x)
        return s if s.startswith("chr") else f"chr{s}"

    def safe_float_mean(cell):
        if pd.isna(cell):
            return np.nan
        vals = []
        for p in re.split(r"[;,]", str(cell).strip()):
            p = p.strip()
            if p and p != "NA":
                try:
                    vals.append(float(p))
                except ValueError:
                    pass
        return float(np.mean(vals)) if vals else np.nan

    def load_mapping(_map_path, _species):
        if _map_path is None:
            file = (pkg_resources.files(mouse_ref).joinpath("mmu_mapping_all.txt")
                    if str(_species).lower() == "mouse"
                    else pkg_resources.files(reference).joinpath("hsa_mapping_all.txt"))
            _map_path = str(file)
        mdf = pd.read_csv(_map_path, sep="\t", usecols=["symbol", "ensembl"])
        mdf["ensembl_norm"] = mdf["ensembl"].map(norm_id)
        return dict(zip(mdf["ensembl_norm"], mdf["symbol"]))

    mapping = load_mapping(map_path, species)

    def to_symbol(gid: str) -> str:
        gidn = norm_id(gid)
        return mapping.get(gidn, gidn)

    def read_rmats(path: str, single: bool = False, inclevel: int = 1) -> pd.DataFrame:
        if single:
            if not os.path.isfile(path):
                raise FileNotFoundError(f"Single rMATS file not found: {path}")
            files = [path]
        else:
            files = sorted(p for p in glob.glob(os.path.join(path, "*")) if os.path.isfile(p))
            if not files:
                raise FileNotFoundError(f"No files found in directory: {path}")
        if inclevel == 1:
            inclevel_str = "IncLevel1"
        else:
            inclevel_str = "IncLevel2"
        def pick_cols(rdf):
            needed = {
                "GeneID": ["GeneID"],
                "geneSymbol": ["geneSymbol"],
                "chr": ["chr"],
                "strand": ["strand"],
                "exonStart": ["exonStart_0base", "exonStart"],
                "exonEnd": ["exonEnd"],
                "upstreamEE": ["upstreamEE"],
                "downstreamES": ["downstreamES"],
                inclevel_str: [inclevel_str, "IncLevel"],
            }
            colmap = {}
            for k, alts in needed.items():
                colmap[k] = next((a for a in alts if a in rdf.columns), None)
            return colmap

        tables = []
        for f in files:
            sample = os.path.splitext(os.path.basename(f))[0]
            rdf = pd.read_csv(f, sep="\t", dtype=str)
            cm = pick_cols(rdf)

            for k, col in cm.items():
                if col is None:
                    rdf[k] = np.nan
                    cm[k] = k

            tmp = pd.DataFrame({
                "GeneID": rdf[cm["GeneID"]].astype(str).map(norm_id),
                "chr": rdf[cm["chr"]],
                "strand": rdf[cm["strand"]],
                "exonStart": pd.to_numeric(rdf[cm["exonStart"]], errors="coerce"),
                "exonEnd": pd.to_numeric(rdf[cm["exonEnd"]], errors="coerce"),
                "upstreamEE": pd.to_numeric(rdf[cm["upstreamEE"]], errors="coerce"),
                "downstreamES": pd.to_numeric(rdf[cm["downstreamES"]], errors="coerce"),
            })
            tmp[sample] = rdf[cm[inclevel_str]].apply(safe_float_mean)

            if cm["geneSymbol"] in rdf.columns:
                gs = rdf[cm["geneSymbol"]].fillna("")
                tmp["geneSymbol"] = [
                    g if (isinstance(g, str) and g.strip()) else to_symbol(ge)
                    for g, ge in zip(gs, tmp["GeneID"])
                ]
            else:
                tmp["geneSymbol"] = tmp["GeneID"].map(to_symbol)

            tables.append(tmp)

        key = ["GeneID", "geneSymbol", "chr", "strand", "exonStart", "exonEnd", "upstreamEE", "downstreamES"]
        out = None
        for t in tables:
            out = t if out is None else pd.merge(out, t, on=key, how="outer")
        return out[key + [c for c in out.columns if c not in key]]

    if single_rMATS:
        df = read_rmats(in_file, single=True)
    else:
        df = read_rmats(in_file, single=False)

    df["chr"] = df["chr"].map(lambda x: fix_chr(x) if pd.notna(x) else x)
    return df


def preprocess_pooled(
    compare_path: str,
    out_psi_dir: str = None,
    *,
    background: Optional[str] = None,
    background_download_root: Optional[str] = None,
    background_path: Optional[str] = None,
    map_path: Optional[str] = None,
    species: str = "human",
    single_rMATS_compare: bool = False,
    single_rMATS_background: bool = False,
    inclevel: int = 1,
) -> Dict[str, Any]:
    """
    Either download a GTEx splicing matrix for `tissue` or use your own `background_path`,
    then combine spliced exons and compute delta PSI.

    Parameters:
        compare_path:   Path to your compare.txt output file. You can also pass in a single rMATS output txt file (SE.MAT.JC.txt or SE.MATS.JCEC.txt), or a folder with rMATS SE.MATS.JC.txt or SE.MATS.JCEC.txt files and the function will process them for you. For rMATS, please rename the files in the folder to the names of your samples.
        out_psi_dir:        Output directory containing PSI files that can be inputed to splitpea.

        tissue:         (optional) Tissue name from IRIS that will be dowloaded. 
                        Y. Pan, J.W. Phillips, B.D. Zhang, M. Noguchi, E. Kutschera, J. McLaughlin, P.A. Nesterenko, Z. Mao, N.J. Bangayan, R. Wang, W. Tran, H.T. Yang, Y. Wang, Y. Xu, M.B. Obusan, D. Cheng, A.H. Lee, K.E. Kadash-Edmondson, A. Champhekar, C. Puig-Saus, A. Ribas, R.M. Prins, C.S. Seet, G.M. Crooks, O.N. Witte, & Y. Xing, IRIS: Discovery of cancer immunotherapy targets arising from pre-mRNA alternative splicing, Proc. Natl. Acad. Sci. U.S.A. 120 (21) e2221116120, https://doi.org/10.1073/pnas.2221116120 (2023).

        background_download_root:  (required if `tissue` is set) Root directory under which to
                        create GTEx_<Tissue>/splicing_matrix/.

        background_path:    (optional) Path to a pre-downloaded background splicing matrix file. You can also pass in a .psi file or a folder with rMATS SE.MATS.JC.txt or SE.MATS.JCEC.txt files and the function will process them for you. For rMATS, please rename the files in the folder to the names of your samples. 
        single_rMATS_compare: (optional) If True, treat the compare file as a single rMATS file.
        single_rMATS_background: (optional) If True, treat the background as a single rMATS file.
        inclevel: (optional) Which rMATS inclusion-level field to use when parsing a single file; 1 or 2. Defaults to 1 if not specified.

    e.g:
        preprocess_pooled(
            tissue="uterus",
            download_root="/uterus",
            compare_path="splicing_matrix.SE.cov10.TCGA_UCS_T.txt",
            psi_dir="/UCS_uterus"
        )
    """

    if out_psi_dir is None:
        out_psi_dir = "out_psi/"

    if bool(background) == bool(background_path):
        raise ValueError("You must provide exactly one of `background` or `background_path`.")
    
    if background_path:
        if not os.path.isfile(background_path):
            raise FileNotFoundError(f"background_path not found: {background_path}")
        
        if os.path.isdir(background_path) or single_rMATS_background: 
            df_bg = read_sample_specific(
                in_file=background_path,
                map_path=map_path,
                species=species,
                single_rMATS=single_rMATS_background,
                inclevel=inclevel
            )
            work_dir = os.path.dirname(background_path) if os.path.isfile(background_path) else background_path
            background_file = os.path.join(work_dir, "background_SE_data.txt")
            df_bg.to_csv(background_file, sep="\t", index=False)
        else:
            background_file = background_path
            work_dir = os.path.dirname(background_path)

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

        slug = background.lower().replace(" ", "").replace("_", "")
        if slug not in slug_map:
            raise ValueError(
                f"Unknown tissue '{background}'. Valid choices are: {', '.join(_ALLOWED_TISSUES)}"
            )
        suffix = slug_map[slug]

        if not background_download_root:
            background_download_root = "IRIS_background/"
            #raise ValueError("`background_download_root` must be set when using bundled `background`.")

        remote_base = (
            "https://xinglabtrackhub.research.chop.edu"
            "/iris/IRIS_data.v2.0.0/db"
        )
        fname       = f"splicing_matrix.SE.cov10.GTEx_{suffix}.txt"
        url         = f"{remote_base}/GTEx_{suffix}/splicing_matrix/{fname}"

        work_dir    = os.path.join(background_download_root, f"GTEx_{suffix}", "splicing_matrix")
        os.makedirs(work_dir, exist_ok=True)
        background_file = os.path.join(work_dir, fname)

        if not os.path.exists(background_file):
            print(f"Downloading {url} to {background_file}")
            print("Y. Pan, J.W. Phillips, B.D. Zhang, M. Noguchi, E. Kutschera, J. McLaughlin, P.A. Nesterenko, Z. Mao, N.J. Bangayan, R. Wang, W. Tran, H.T. Yang, Y. Wang, Y. Xu, M.B. Obusan, D. Cheng, A.H. Lee, K.E. Kadash-Edmondson, A. Champhekar, C. Puig-Saus, A. Ribas, R.M. Prins, C.S. Seet, G.M. Crooks, O.N. Witte, & Y. Xing, IRIS: Discovery of cancer immunotherapy targets arising from pre-mRNA alternative splicing, Proc. Natl. Acad. Sci. U.S.A. 120 (21) e2221116120, https://doi.org/10.1073/pnas.2221116120 (2023).")
            resp = requests.get(url, stream=True)
            resp.raise_for_status()
            with open(background_file, "wb") as wf:
                for chunk in resp.iter_content(chunk_size=8192):
                    wf.write(chunk)

    combine_spliced_exon(work_dir)

    combined_list = glob.glob(os.path.join(work_dir, "*_combined_mean.txt"))
    if not combined_list:
        raise FileNotFoundError(f"No *_combined_mean.txt found in {work_dir}")
    
    base = os.path.splitext(os.path.basename(background_file))[0]
    combined_file = os.path.join(work_dir, base + "_combined_mean.txt")
    if not os.path.exists(combined_file):
        print(
            f"Warning: expected combined file '{combined_file}' was not found. "
            f"Falling back to '{combined_list[0]}'."
        )
        combined_file = combined_list[0]

    if os.path.isdir(compare_path) or single_rMATS_compare:
        df_cmp = read_sample_specific(
            in_file=compare_path,
            map_path=map_path,
            species=species,
            single_rMATS=single_rMATS_compare,
            inclevel=inclevel
        )
        work_dir_ = os.path.dirname(compare_path) if os.path.isfile(compare_path) else compare_path
        compare_path = os.path.join(work_dir_, "compare_SE_data.txt")
        df_cmp.to_csv(compare_path, sep="\t", index=False)

    print("Starting delta PSI calculation...")
    delta = calculate_delta_psi(
        combined_file,  
        background_file,    
        compare_path,
        out_psi_dir
    )

    return out_psi_dir