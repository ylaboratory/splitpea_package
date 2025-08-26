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

    # 1) Load and preprocess background_samples summary (one PSI column)
    sum_bg = pd.read_csv(sum_bg_file, sep='\t')
    sum_bg = sum_bg.groupby(group_cols, as_index=False).mean()
    psi_cols = [c for c in sum_bg.columns if c not in group_cols]
    if len(psi_cols) != 1:
        raise ValueError("Expected exactly one PSI column in summarized background file.")
    sum_bg = sum_bg.rename(columns={psi_cols[0]: 'psi_background_samples'})

    # 2) Load and preprocess raw background_samples for ECDF
    bg_all = pd.read_csv(bg_file, sep='\t')
    bg_all = bg_all.rename(columns={'exonStart':'exonStart_0base'})
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

    # 3) Load and preprocess compare_samples data
    compare_samples_all = pd.read_csv(target_file, sep='\t')
    compare_samples_all = compare_samples_all.rename(columns={'GeneID':'AC','geneSymbol':'GeneName'})
    compare_samples_all = compare_samples_all.drop(columns=['upstreamEE','downstreamES'], errors='ignore')
    compare_samples_all = compare_samples_all.groupby(['AC','GeneName','chr','strand','exonStart','exonEnd'], as_index=False).mean()
    compare_samples_all['chr'] = compare_samples_all['chr'].astype(str).str.replace('chr','')
    compare_samples_all['strand'] = compare_samples_all['strand'].astype(str).map(lambda s: -1 if s=='-' else 1)
    compare_samples_all = compare_samples_all.rename(columns={
        'AC':'ensembl_gene_id',
        'GeneName':'symbol',
        'exonStart':'exon_start',
        'exonEnd':'exon_end'
    })
    sample_cols_compare_samples = [c for c in compare_samples_all.columns if c not in group_cols]
    compare_samples_all = compare_samples_all[group_cols + sample_cols_compare_samples]

    # 4) For each compare_samples sample: compute delta PSI and p-values
    for sam in sample_cols_compare_samples:
        print(f"Processing sample {sam}...")
        tmp = compare_samples_all[group_cols + [sam]].copy()
        tmp = tmp.rename(columns={sam: 'psi_compare_samples'})
        
        # merge with background_samples summary
        joint = pd.merge(sum_bg, tmp, on=group_cols, how='inner')
        joint = joint.dropna(subset=['psi_background_samples','psi_compare_samples'])
        joint['delta_psi'] = joint['psi_compare_samples'] - joint['psi_background_samples']
        
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


def read_sample_specific(
    in_file: str,
    map_path: str = None,
    species: str = "Human",
    gtf: str = None,
) -> pd.DataFrame:
    """
    Convert SUPPA2 .psi file(s) or rMATS output into a unified PSI table.

    Parameters:
    - in_file: path to a SUPPA2 .psi file or a directory of rMATS files
    - map_path: optional Ensembl→symbol mapping TSV with columns [symbol, ensembl]
    - species: "Human" or "mouse" (used if map_path not given)
    - gtf: optional GTF file (only used for SUPPA2) to recover skipped exon start/end

    Outputs:
    - DataFrame columns:
      [GeneID, geneSymbol, chr, strand, exonStart, exonEnd,
       upstreamEE, downstreamES, sample1, sample2, ...]
    """

    # internal helper functions
    def norm_id(x: str) -> str:
        if pd.isna(x): return x
        s = str(x)
        return s.split(".", 1)[0].split("_", 1)[0]

    def fix_chr(x: str):
        if pd.isna(x): return x
        s = str(x)
        return s if s.startswith("chr") else f"chr{s}"

    def safe_float_mean(cell):
        if pd.isna(cell): return np.nan
        vals = []
        for p in re.split(r"[;,]", str(cell).strip()):
            p = p.strip()
            if p and p != "NA":
                try: vals.append(float(p))
                except ValueError: pass
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

    def load_gtf_index(gtf_path: str):
        if gtf_path is None or not os.path.isfile(gtf_path):
            return None
        cols = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
        g = pd.read_csv(
            gtf_path, sep="\t", comment="#", header=None, names=cols,
            dtype={"seqname": str, "feature": str, "start": int, "end": int, "strand": str, "attribute": str}
        )
        g = g[g["feature"] == "exon"].copy()

        def get_attr(attr, key):
            m = re.search(rf'{key}\s+"([^"]+)"', attr)
            return m.group(1) if m else None

        g["gene_id"] = g["attribute"].apply(lambda a: norm_id(get_attr(a, "gene_id")))
        g["seqname"] = g["seqname"].map(fix_chr)
        idx = {}
        for (gid, ch, st), blk in g.groupby(["gene_id", "seqname", "strand"]):
            idx[(gid, ch, st)] = blk[["start", "end"]].sort_values(["start", "end"]).to_numpy()
        return idx

    gtf_index = load_gtf_index(gtf)

    def fallback_exon(a2, b1):
        if pd.notna(a2) and pd.notna(b1) and (b1 - a2 >= 2):
            return int(a2 + 1), int(b1 - 1)
        return np.nan, np.nan

    def infer_exon_from_gtf(gid, ch, st, a2, b1):
        if gtf_index is None or any(pd.isna(v) for v in [gid, ch, st, a2, b1]):
            return fallback_exon(a2, b1)
        exons = gtf_index.get((str(gid), str(ch), str(st)))
        if exons is None: return fallback_exon(a2, b1)
        L, R = int(a2), int(b1)
        best, best_len = None, -1
        for start, end in exons:
            if start > L and end < R:
                length = end - start + 1
                if length > best_len:
                    best, best_len = (int(start), int(end)), length
        if best is not None:
            return best
        for start, end in exons:
            if start >= L and end <= R:
                return int(start), int(end)
        return fallback_exon(a2, b1)

    def parse_suppa_event(event_id: str):
      """
      SE-only parser.
      Example: ENSG...;SE:chr11:125496728-125497502:125497725-125499127:+
      """
      s = str(event_id)

      m = re.search(r";SE:[Cc][Hh][Rr]?([^:]*):(\d+)-(\d+):(\d+)-(\d+):([+-])", s)
      if not m:
          # Be tolerant in case some tools drop the leading semicolon:
          m = re.search(r"\bSE:[Cc][Hh][Rr]?([^:]*):(\d+)-(\d+):(\d+)-(\d+):([+-])", s)
          if not m:
              return None  

      gid = norm_id(s.split(";", 1)[0])

      chr_part = m.group(1)
      a1, a2 = int(m.group(2)), int(m.group(3))  # upstream exon start-end
      b1, b2 = int(m.group(4)), int(m.group(5))  # downstream exon start-end
      strand = m.group(6)

      return {
          "GeneID": gid,
          "chr": fix_chr(chr_part),
          "strand": strand,
          "a2": a2, "b1": b1,
          "upstreamEE": a2,
          "downstreamES": b1,
      }

    def read_suppa(path: str) -> pd.DataFrame:
      import io

      rows = []
      with open(path, "r", encoding="utf-8") as fh:
          for line in fh:
              line = line.strip()
              if not line:
                  continue
              parts = re.split(r"\s+", line)
              rows.append(parts)

      if not rows:
          raise ValueError("Empty SUPPA2 file.")

      def looks_like_event(s: str) -> bool:
          return isinstance(s, str) and (";SE:" in s or s.startswith(("ENSG", "ENSMUSG")))

      first_cell = rows[0][0]
      if looks_like_event(first_cell):
          ncols = len(rows[0])
          cols = ["Event_ID"] + [f"sample{i}" for i in range(1, ncols)]
          data = pd.DataFrame(rows, columns=cols)
      else:
          hdr = rows[0]
          data = pd.DataFrame(rows[1:])
          if len(hdr) == data.shape[1] - 1:
              hdr = ["Event_ID"] + hdr
          if len(hdr) != data.shape[1]:
              hdr = ["Event_ID"] + [f"sample{i}" for i in range(1, data.shape[1])]
          data.columns = hdr
          if data.columns[0].lower() not in ("event_id", "event", "id"):
              data.rename(columns={data.columns[0]: "Event_ID"}, inplace=True)

      if "Event_ID" not in data.columns:
          data = data.reset_index().rename(columns={"index": "Event_ID"})

      sample_cols = [c for c in data.columns if c != "Event_ID"]
      meta_series = data["Event_ID"].apply(parse_suppa_event)
      keep_mask = meta_series.notnull()
      if not keep_mask.any():
          raise ValueError("No SE events found in SUPPA2 file.")

      meta = meta_series[keep_mask].apply(pd.Series).reset_index(drop=True)
      samples = data.loc[keep_mask, sample_cols].apply(pd.to_numeric, errors="coerce").reset_index(drop=True)

      meta["geneSymbol"] = meta["GeneID"].apply(to_symbol)

      if gtf_index is not None:
          exon_bounds = meta.apply(
              lambda r: infer_exon_from_gtf(r["GeneID"], r["chr"], r["strand"], r["a2"], r["b1"]),
              axis=1, result_type="expand"
          )
      else:
          exon_bounds = meta.apply(
              lambda r: fallback_exon(r["a2"], r["b1"]),
              axis=1, result_type="expand"
          )
      exon_bounds.columns = ["exonStart", "exonEnd"]

      out = pd.concat(
          [meta[["GeneID", "geneSymbol", "chr", "strand", "upstreamEE", "downstreamES"]],
          exon_bounds, samples],
          axis=1
      )

      ordered = ["GeneID", "geneSymbol", "chr", "strand",
                "exonStart", "exonEnd", "upstreamEE", "downstreamES"] + sample_cols
      return out[ordered]

    def read_rmats_dir(path: str) -> pd.DataFrame:
        files = sorted(p for p in glob.glob(os.path.join(path, "*")) if os.path.isfile(p))
        if not files:
            raise FileNotFoundError(f"No files found in directory: {path}")

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
                "IncLevel1": ["IncLevel1", "IncLevel"],
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
            tmp[sample] = rdf[cm["IncLevel1"]].apply(safe_float_mean)

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

    df = read_rmats_dir(in_file) if os.path.isdir(in_file) else read_suppa(in_file)
    df["chr"] = df["chr"].map(lambda x: fix_chr(x) if pd.notna(x) else x)
    return df


def preprocess_pooled(
    compare_file: str,
    out_psi_dir: str,
    *,
    background: Optional[str] = None,
    background_download_root: Optional[str] = None,
    background_path: Optional[str] = None,
    map_path: Optional[str] = None,
    species: str = "Human",
    gtf: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Either download a GTEx splicing matrix for `tissue` or use your own `background_path`,
    then combine spliced exons and compute delta PSI.

    Args:
        compare_file:   Path to your compare.txt output file. You can also pass in a .psi file or a folder with rMATS SE.MATS.JC.txt or SE.MATS.JCEC.txt files and the function will process them for you. For rMATS, please rename the files in the folder to the names of your samples. 
        out_psi_dir:        Output directory containing PSI files that can be inputed to splitpea.

        tissue:         (optional) Tissue name from IRIS that will be dowloaded. 
                        Y. Pan, J.W. Phillips, B.D. Zhang, M. Noguchi, E. Kutschera, J. McLaughlin, P.A. Nesterenko, Z. Mao, N.J. Bangayan, R. Wang, W. Tran, H.T. Yang, Y. Wang, Y. Xu, M.B. Obusan, D. Cheng, A.H. Lee, K.E. Kadash-Edmondson, A. Champhekar, C. Puig-Saus, A. Ribas, R.M. Prins, C.S. Seet, G.M. Crooks, O.N. Witte, & Y. Xing, IRIS: Discovery of cancer immunotherapy targets arising from pre-mRNA alternative splicing, Proc. Natl. Acad. Sci. U.S.A. 120 (21) e2221116120, https://doi.org/10.1073/pnas.2221116120 (2023).

        background_download_root:  (required if `tissue` is set) Root directory under which to
                        create GTEx_<Tissue>/splicing_matrix/.

        background_path:    (optional) Path to a pre-downloaded background splicing matrix file. You can also pass in a .psi file or a folder with rMATS SE.MATS.JC.txt or SE.MATS.JCEC.txt files and the function will process them for you. For rMATS, please rename the files in the folder to the names of your samples. 
    
    e.g:
        preprocess_pooled(
            tissue="uterus",
            download_root="/uterus",
            compare_file="splicing_matrix.SE.cov10.TCGA_UCS_T.txt",
            psi_dir="/UCS_uterus"
        )
    """
    if bool(background) == bool(background_path):
        raise ValueError("You must provide exactly one of `background` or `background_path`.")
    
    if background_path:
        if not os.path.isfile(background_path):
            raise FileNotFoundError(f"background_path not found: {background_path}")
        
        if background_path.endswith(".psi") or os.path.isdir(background_path): # if it is SUPPA2 (.psi) or a folder of rMATS files
            df_bg = read_sample_specific(
                in_file=background_path,
                map_path=map_path,
                species=species,
                gtf=gtf,
            )
            work_dir = os.path.dirname(background_path) if os.path.isfile(background_path) else background_path
            background_file = os.path.join(work_dir, "background.psi.txt")
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
            raise ValueError("`background_download_root` must be set when using bundled `background`.")

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

    if compare_file.endswith(".psi") or os.path.isdir(compare_file):
        df_cmp = read_sample_specific(
            in_file=compare_file,
            map_path=map_path,
            species=species,
            gtf=gtf,
        )
        os.makedirs(out_psi_dir, exist_ok=True)
        compare_file = os.path.join(out_psi_dir, "compare.psi.txt")
        df_cmp.to_csv(compare_file, sep="\t", index=False)

    print("Starting delta PSI calculation...")
    delta = calculate_delta_psi(
        combined_file,  
        background_file,    
        compare_file,
        out_psi_dir
    )

    return out_psi_dir