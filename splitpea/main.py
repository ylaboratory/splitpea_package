# splitpea/main.py

import argparse
import pickle
from collections import defaultdict
import networkx as nx
import os, sys
import numpy as np
import pandas as pd
from statistics import mean
from functools import partial
from os.path import basename
from subprocess import Popen, PIPE
from pathlib import Path
import logging
import csv
import tempfile
import shutil

from .exons import Exons
from .parser import parse_suppa2
from .parser import parse_rmats
from .grab_stats import rewired_edges_stat
from .grab_stats import rewired_genes
from .get_background_ppi import get_background
from .plot_network import plot_rewired_network
from .core import tb_query
from .core import calculate_edges
from.core import deduce_final_edge_weights

import importlib_resources as pkg_resources
from .src import reference
from .src import mouse_ref

logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s',
                    datefmt="%m-%d-%Y %I:%M:%S %p",
                    level=logging.ERROR)
logger = logging.getLogger('diff_exon')


def calculate_delta_psi(sum_bg_file, bg_file, target_file, outdir):
    """
    Calculate delta PSI values and p-values for GTEx–TCGA exon comparisons.
    
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

def fileExists(f):
    """
    Check if a file exists.
    """
    if not os.path.isfile(f):
        print("ERROR: the file " + f + " was not found.")


def run(in_file,
        out_file_prefix: str,
        skip: int = 1,
        dpsi_cut: float = 0.05,
        sigscore_cut: float = 0.05,
        include_nas: bool = True,
        verbose: bool = False,
        input_format: str = "regular",
        ppif: str = None,
        ddif: str = None,
        entrezpfamf: str = None,
        pfamcoordsf: str = None,
        tbf: str = None,
        species: str = 'human',
        index: int = None,
        edge_stats_file: str = None,
        gene_degree_stats: bool = False,
        plot_net: bool = False,
        gephi_tsv: bool = False,
        cytoscape_gml: bool = False,
        map_path: str = None):
    """
    Run the splicing-specific network pipeline.

    Parameters:
      in_file: If input_format is 'regular', a file path (str) to the input file with
               differentially expressed exons.
               If input_format is 'suppa2', a list/tuple of 2 file paths: [psivec_file, dpsi_file].
      out_file_prefix: Prefix for output files.
      skip: Number of lines to skip in the input file (default=0).
      dpsi_cut: Delta PSI cutoff (default=0.05).
      sigscore_cut: Significance score cutoff (default=0.05).
      include_nas: Include NAs in significance testing (default=True).
      verbose: Enable verbose logging (default=False).
      input_format: Either 'regular' (default) or 'suppa2'.
      species: Defaults to human.
      index: If exon files has indx 1 or 0-based. 
      edge_stats_file: Path to file that stores rewired edge stats. 
      gene_degree_stats: Path to file to store rewired gene degree stats.
      plot_net: Plots rewired network (default = False).
      
      The following file paths default to files under 'src/reference':
      ppif: Protein-protein interaction file (default: human_ppi_0.5.dat).
      ddif: Domain-domain interaction file (default: ddi_0.5.dat).
      entrezpfamf: Gene-protein domain info file (default: human_entrez_pfam.txt).
      pfamcoordsf: Pfam genome coordinates file (default: human_pfam_genome_coords_sorted.txt.gz).
      tbf: Tabix file (default: human_pfam_genome_coords_sorted.txt.gz).
      map_path: Mapping between gene ids (default: hsa_mapping_all.txt). 
      Note: there are equivalents in 'src/mouse_ref'
    """
    if verbose:
        logger.setLevel(logging.INFO)

    if input_format.lower() == "suppa2":
        if not (isinstance(in_file, (list, tuple)) and len(in_file) == 2):
            raise ValueError("For SUPPA2 input format, in_file must be a list or tuple of 2 file paths: [psivec_file, dpsi_file].")
        psivec_file, dpsi_file = in_file
        
        if psivec_file.endswith(".dpsi") and dpsi_file.endswith(".psivec"):
            psivec_file, dpsi_file = dpsi_file, psivec_file  # Swap variables

        if map_path is None:
            if species == "mouse":
                file = pkg_resources.files(mouse_ref).joinpath("mmu_mapping_all.txt")
            else:
                file = pkg_resources.files(reference).joinpath("hsa_mapping_all.txt")
            map_path = str(file)  
        fileExists(map_path)
        in_file = parse_suppa2(psivec_file, dpsi_file, map_path, species=species, verbose = verbose)

        if skip is None:
            skip = 1
        if index is None:
            index = 1
    else:
        if type(in_file) != str:
            if len(in_file) > 1:
                raise ValueError("For input_formats other than SUPPA2, only a single input file is allowed.")
            in_file = in_file[0]
        if index is None:
            index = 0 
    
    if input_format.lower() == "rmats":
        in_file = parse_rmats(in_file)
        if skip is None:
            skip = 1

    # Set default reference file paths if not provided.
    if ppif is None or ddif is None or entrezpfamf is None or pfamcoordsf is None or tbf is None:
        if species.lower() == "mouse":
            if ppif is None:
                ppif = str(pkg_resources.files(mouse_ref).joinpath("mouse_ppi.dat"))
            if ddif is None:
                ddif = str(pkg_resources.files(mouse_ref).joinpath("ddi_0.5.dat"))
            if entrezpfamf is None:
                entrezpfamf = str(pkg_resources.files(mouse_ref).joinpath("mouse_entrez_pfam.txt"))
            if pfamcoordsf is None:
                pfamcoordsf = str(pkg_resources.files(mouse_ref).joinpath("mouse_pfam_genome_coords_sorted.txt.gz"))
            if tbf is None:
                tbf = str(pkg_resources.files(mouse_ref).joinpath("mouse_pfam_genome_coords_sorted.txt.gz"))
        else:
            if ppif is None:
                ppif = str(pkg_resources.files(reference).joinpath("human_ppi_0.5.dat"))
            if ddif is None:
                ddif = str(pkg_resources.files(reference).joinpath("ddi_0.5.dat"))
            if entrezpfamf is None:
                entrezpfamf = str(pkg_resources.files(reference).joinpath("human_entrez_pfam.txt"))
            if pfamcoordsf is None:
                pfamcoordsf = str(pkg_resources.files(reference).joinpath("human_pfam_genome_coords_sorted.txt.gz"))
            if tbf is None:
                tbf = str(pkg_resources.files(reference).joinpath("human_pfam_genome_coords_sorted.txt.gz"))

    # Check that reference files exist.
    fileExists(ppif)
    fileExists(entrezpfamf)
    fileExists(pfamcoordsf)
    fileExists(pfamcoordsf + ".tbi")

    logger.info("Loading DDIs....")
    all_ddis = nx.read_weighted_edgelist(ddif)
    logger.info("# pfams: %i, # ddis: %i", len(all_ddis), all_ddis.number_of_edges())

    logger.info("Loading PPIs....")
    all_ppis = nx.read_weighted_edgelist(ppif)
    logger.info("# proteins: %i, # interactions: %i", len(all_ppis), all_ppis.number_of_edges())

    logger.info("Loading gene-protein domain info....")
    entrez_pfams = nx.bipartite.read_edgelist(entrezpfamf)

    logger.info("Reading differential exon results....")
    dexons = Exons()
    dexons.read_dexons(in_file, index, skip, dpsi_cut, sigscore_cut, include_nas)

    logger.info("# of differentially expressed exons: %i", len(dexons))

    logger.info("Calculating network....")
    diff_splice_g = nx.Graph(name=basename(in_file))
    calculate_edges(diff_splice_g, dexons, tbf, all_ddis, all_ppis, entrez_pfams)
    deduce_final_edge_weights(diff_splice_g)

    logger.info("# nodes: %i, # edges: %i", len(diff_splice_g), diff_splice_g.number_of_edges())

    logger.info("Outputting as pickle...")
    with open(out_file_prefix + '.edges.pickle', 'wb') as out:
        pickle.dump(diff_splice_g, out)

    with open(out_file_prefix + '.edges.dat', 'w') as out:
        header = ["node1", "node2", "weight", "chaos"]
        out.write('\t'.join(header) + '\n')
        for gi, gj in nx.edges(diff_splice_g):
            ppi_ddi_out = [gi, gj, str(diff_splice_g[gi][gj]['weight']),
                           str(diff_splice_g[gi][gj]['chaos'])]
            out.write('\t'.join(ppi_ddi_out) + '\n')
    
    if edge_stats_file is not None:
        dat_file = out_file_prefix + '.edges.dat'

        stats = rewired_edges_stat(dat_file)

        file_exists = os.path.exists(edge_stats_file)
        with open(edge_stats_file, 'a') as f:
            if not file_exists:
                f.write("sample\tnum_gain_edges\tnum_loss_edges\tnum_chaos_edges\n")
                logger.info("Edge stat file created: " + str(edge_stats_file))
            f.write(f"{out_file_prefix}\t{stats['gain']}\t{stats['loss']}\t{stats['chaos']}\n")
            logger.info("Edge stat file written")

    if gene_degree_stats == True:
        background = get_background(ppif, ddif, entrezpfamf)
        gene_stats_df = rewired_genes(diff_splice_g, background, map_path)
        gene_stats_df.to_csv(out_file_prefix + "_gene_degree.csv", index=False)

    if plot_net == True:
        logger.info("Plotting network...")
        plot_rewired_network(diff_splice_g, with_labels = True, pdf_path=out_file_prefix + "_network_plot.pdf")
    if gephi_tsv == True:
        logger.info("Outputing gexf for network...")
        plot_rewired_network(diff_splice_g, with_labels = True, plot_matplot = False, gephi_path=out_file_prefix + '_gephi.csv')
    if cytoscape_gml == True:
        logger.info("Outputing gml for network...")
        plot_rewired_network(diff_splice_g, with_labels = True, plot_matplot = False, cytoscape_path=out_file_prefix + '_cyto.gml')
    
    if input_format.lower() in ("rmats", "suppa2"):
        try:
            os.remove(in_file)
            logger.info("Temporary file %s deleted.", in_file)
        except Exception as e:
            logger.error("Failed to delete temporary file %s: %s", in_file, e)

    logger.info("Done")

    return diff_splice_g

def plot(pickle_path: str,
         with_labels: bool = False,
         pdf_path: str = None,
         gephi_path: str = None,
         cytoscape_path: str = None,
         max_nodes: int = 2000,
         max_edges: int = 10000,
         symbol: bool = True,
         map_path: str = None,
         species: str = 'human'):
    """
    Load a pickled Graph (as created by `rewire(...)`) and call plot_rewired_network().

    Parameters:
    - pickle_path: path to the '.edges.pickle' file (output of `rewire(...)`).
    - with_labels: whether to draw node labels in the plot.
    - pdf_path: if provided, write a PDF of the plotted network to this path.
    - gephi_path: if provided, write a Gephi-compatible CSV to this path.
    - cytoscape_path: if provided, write a Cytoscape GML to this path.
    """

    if not os.path.isfile(pickle_path):
        sys.exit(f"ERROR: pickle file '{pickle_path}' not found.")

    with open(pickle_path, 'rb') as f:
        G = pickle.load(f)
    
    if symbol:
        if map_path is None or not os.path.isfile(map_path):
            if species == "mouse":
                file = pkg_resources.files(mouse_ref).joinpath("mmu_mapping_all.txt")
            else:
                file = pkg_resources.files(reference).joinpath("hsa_mapping_all.txt")
            map_path = str(file)  

        mapping_df = pd.read_csv(map_path, sep="\t", dtype=str)
        entrez_map  = mapping_df.set_index('entrez')['symbol'].to_dict()
        ensembl_map = mapping_df.set_index('ensembl')['symbol'].to_dict()
        uniprot_map = mapping_df.set_index('uniprot')['symbol'].to_dict()

        relabel_map = {}
        for node in G.nodes():
            s = None
            node_str = str(node)
            s = entrez_map.get(node_str) or ensembl_map.get(node_str) or uniprot_map.get(node_str)
            relabel_map[node] = s if s is not None else node  

        G = nx.relabel_nodes(G, relabel_map)

    plot_matplot = True
    if pdf_path is None:
        plot_matplot = False

    result = plot_rewired_network(
        G,
        with_labels=with_labels,
        pdf_path=pdf_path,
        plot_matplot=plot_matplot,
        gephi_path=gephi_path,
        cytoscape_path=cytoscape_path,
        max_nodes = max_nodes,
        max_edges = max_edges
    )

    outputs = []
    if pdf_path:
        outputs.append(f"PDF → {pdf_path}")
    if gephi_path:
        outputs.append(f"Gephi TSV → {gephi_path}")
    if cytoscape_path:
        outputs.append(f"Cytoscape GML → {cytoscape_path}")

    if outputs:
        print("Generated outputs:\n  " + "\n  ".join(outputs))
    else:
        print("No output files were requested.")

def stats(dat_file: str,
          rewire_net: str,
          out_file_prefix: str = None,
          ppif: str = None,
          ddif: str = None,
          entrezpfamf: str = None,
          map_path: str = None, 
          species: str = 'human'):
    """
    Output gene-level statistics for a rewired splicing-induced PPI network.

    Parameters:
      dat_file: Path to the input file containing the rewired network edges (e.g., one edge per line).
      rewire_net: Path to the file listing genes with differential splicing events.
      out_file: Path (or prefix) for the output gene‐stats summary.
      ppif: Protein–protein interaction reference file (default: human_ppi_0.5.dat).
      ddif: Domain–domain interaction reference file (default: ddi_0.5.dat).
      entrezpfamf: File mapping Entrez Gene IDs to Pfam domains (default: human_entrez_pfam.txt).
      species: Species identifier for reference selection (default: 'human').
    """
    
    if ppif is None or ddif is None or entrezpfamf is None or map_path is None:
        if species.lower() == "mouse":
            if ppif is None:
                ppif = str(pkg_resources.files(mouse_ref).joinpath("mouse_ppi.dat"))
            if ddif is None:
                ddif = str(pkg_resources.files(mouse_ref).joinpath("ddi_0.5.dat"))
            if entrezpfamf is None:
                entrezpfamf = str(pkg_resources.files(mouse_ref).joinpath("mouse_entrez_pfam.txt"))
            if map_path is None:
                file = pkg_resources.files(mouse_ref).joinpath("mmu_mapping_all.txt")
                map_path = str(file) 
            # if pfamcoordsf is None:
            #     pfamcoordsf = str(pkg_resources.files(mouse_ref).joinpath("mouse_pfam_genome_coords_sorted.txt.gz"))
            # if tbf is None:
            #     tbf = str(pkg_resources.files(mouse_ref).joinpath("mouse_pfam_genome_coords_sorted.txt.gz"))
        else:
            if ppif is None:
                ppif = str(pkg_resources.files(reference).joinpath("human_ppi_0.5.dat"))
            if ddif is None:
                ddif = str(pkg_resources.files(reference).joinpath("ddi_0.5.dat"))
            if entrezpfamf is None:
                entrezpfamf = str(pkg_resources.files(reference).joinpath("human_entrez_pfam.txt"))
            # if pfamcoordsf is None:
            #     pfamcoordsf = str(pkg_resources.files(reference).joinpath("human_pfam_genome_coords_sorted.txt.gz"))
            # if tbf is None:
            #     tbf = str(pkg_resources.files(reference).joinpath("human_pfam_genome_coords_sorted.txt.gz"))
            if map_path is None:
                file = pkg_resources.files(reference).joinpath("hsa_mapping_all.txt")
                map_path = str(file)  


    stats = rewired_edges_stat(dat_file)
    gain, loss, chaos = (int(stats[k]) for k in ("gain","loss","chaos"))
    print(f"""\
    Basic edge stats:
    Gain : {gain}
    Loss : {loss}
    Chaos: {chaos}
    """)

    background = get_background(ppif, ddif, entrezpfamf)
    with open(rewire_net, 'rb') as f:
        diff_splice_g = pickle.load(f)
    gene_stats_df = rewired_genes(diff_splice_g, background, map_path)
    if out_file_prefix != None:
        gene_stats_df.to_csv(out_file_prefix + "_gene_degree.csv", index=False)
    return gene_stats_df


def main():

    parser = argparse.ArgumentParser(
        description="splitpea: SPLicing InTeractions PErsonAlized"
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    rewire_p = subparsers.add_parser(
        "run", help="calculate and output rewired network"
    )

    rewire_p.add_argument(
        'in_file', nargs='+',
        help="Input file path(s). For SUPPA2 format, please provide two files: psivec_path and dpsi_path."
    )
    rewire_p.add_argument('out_file_prefix', 
                        help="Prefix for output files.")
    rewire_p.add_argument('--input_format', type=str, choices=['regular', 'suppa2', 'rmats'], default='regular', 
                        help="Input file format (default: regular).")
    rewire_p.add_argument('--skip', type=int, default=1, 
                        help="Number of lines to skip in the input file (default: 1).")
    rewire_p.add_argument('--dpsi_cut', type=float, default=0.05, 
                        help="Delta PSI cutoff (default: 0.05).")
    rewire_p.add_argument('--sigscore_cut', type=float, default=0.05, 
                        help="Significance score cutoff (default: 0.05).")
    rewire_p.add_argument('--include_nas', type=bool, default=True,  
                        help="Include NAs in significance testing.")
    rewire_p.add_argument('--verbose', action='store_true', 
                        help="Enable verbose logging.")
    rewire_p.add_argument('--ppif', type=str, default=None, 
                        help="Protein-protein interaction file path.")
    rewire_p.add_argument('--ddif', type=str, default=None, 
                        help="Domain-domain interaction file path.")
    rewire_p.add_argument('--entrezpfamf', type=str, default=None, 
                        help="Gene-protein domain info file path.")
    rewire_p.add_argument('--pfamcoordsf', type=str, default=None, 
                        help="Pfam genome coordinates file path.")
    rewire_p.add_argument('--tbf', type=str, default=None, 
                        help="Tabix file path.")
    rewire_p.add_argument('--species', type=str, default='human', 
                        help="Species (default: human).")
    rewire_p.add_argument('--index', type=int, default=None, choices=[0, 1], 
                        help="indexing scheme of data, 0 or 1 (default=0)")
    rewire_p.add_argument('--edge_stats_file', type=str, default=None, 
                    help="File path to a txt file to save Splitpea network statistics, including counts of gained, lost, and chaotic edges. If not file path exists it creates the file.")
    rewire_p.add_argument('--gene_degree_stats', action='store_true', 
                        help="Outputs degree stats of genes in the rewired network and saves a background PPI network.")
    rewire_p.add_argument('--plot_net', action='store_true', 
                        help="Plots a rough version of the rewired network using matplotlib and saves a pdf. For large networks, this may crash. Consider using Gephi, Cytoscape, or another software to plot.")
    rewire_p.add_argument('--gephi_tsv', action='store_true', 
                        help="Saves a .tsv file for more detailed plotting that is compatible with the Gephi software.")
    rewire_p.add_argument('--cytoscape_gml', action='store_true', 
                    help="Saves a .gml file for more detailed plotting that is compatible with the Cytoscape software.")
    rewire_p.add_argument('--map_path', type=str, default=None,
                        help="Path to text file that mapps between gene ids, where it has tab delineated columns: symbol, entrez, ensembl, uniprot. The pacakage has default files for mouse and human.")



    plot_p = subparsers.add_parser(
        "plot", help="load a saved rewired network (.edges.pickle) and produce plots"
    )
    plot_p.add_argument(
        'pickle_path',
        help="Path to the '.edges.pickle' file (output of `rewire`)."
    )
    plot_p.add_argument(
        '--with_labels', action='store_true',
        help="Draw node labels in the plot (default: False)."
    )
    plot_p.add_argument(
        '--pdf_path', type=str, default=None,
        help="If provided, write a PDF plot to this location."
    )
    plot_p.add_argument(
        '--gephi_path', type=str, default=None,
        help="If provided, write a Gephi-compatible CSV to this location."
    )
    plot_p.add_argument(
        '--cytoscape_path', type=str, default=None,
        help="If provided, write a Cytoscape .gml to this location."
    )

    plot_p.add_argument(
        '--max_nodes', type=int, default=2000,
        help="Maximum number of nodes before disabling Matplotlib plotting (default: 2000)."
    )
    plot_p.add_argument(
        '--max_edges', type=int, default=10000,
        help="Maximum number of edges before disabling Matplotlib plotting (default: 10000)."
    )
    plot_p.add_argument(
        '--symbol', type=bool, default=True,
        help="Replace node id with gene symbol"
    )
    plot_p.add_argument(
        '--map_path', type=str, default=None,
        help="Path to text file that mapps between gene ids, where it has tab delineated columns: symbol, entrez, ensembl, uniprot. The pacakage has default files for mouse and human."
    )
    plot_p.add_argument(
        '--species', type=str, default='human',
        help="Species (default: human)."
    )

    stats_p = subparsers.add_parser(
    "stats",
    help="compute and write gene level statistics for a rewired splicing PPI network"
    )
    stats_p.add_argument(
        "dat_file",
        type=str,
        help="Path to the input file of rewired network edges dat file"
    )
    stats_p.add_argument(
        "rewire_net",
        type=str,
        help="Path to pickle file of rewired network" 
    )
    stats_p.add_argument(
        "out_file_prefix",
        type=str,
        help="Prefix for path for the output gene rewire statistics summary"
    )
    stats_p.add_argument(
        "--ppif",
        type=str,
        default=None,
        help="Protein protein interaction reference file (default: from src/reference for the given species)"
    )
    stats_p.add_argument(
        "--ddif",
        type=str,
        default=None,
        help="Domain domain interaction reference file (default: from src/reference for the given species)"
    )
    stats_p.add_argument(
        "--entrezpfamf",
        type=str,
        default=None,
        help="Entrez to Pfam mapping file (default: from src/reference for the given species)"
    )
    stats_p.add_argument(
        "--map_path",
        type=str,
        default=None,
        help="Path to text file that mapps between gene ids, where it has tab delineated columns: symbol, entrez, ensembl, uniprot. The pacakage has default files for mouse and human."
    )
    stats_p.add_argument(
        "--species",
        type=str,
        default="human",
        help="Species identifier to pick the correct reference files (default: 'human')"
    )


    args = parser.parse_args()

    if args.command == "run":
        run(args.in_file,
            args.out_file_prefix,
            skip=args.skip,
            dpsi_cut=args.dpsi_cut,
            sigscore_cut=args.sigscore_cut,
            include_nas=args.include_nas,
            verbose=args.verbose,
            input_format=args.input_format,
            ppif=args.ppif,
            ddif=args.ddif,
            entrezpfamf=args.entrezpfamf,
            pfamcoordsf=args.pfamcoordsf,
            tbf=args.tbf,
            species=args.species,
            index = args.index, 
            edge_stats_file = args.edge_stats_file,
            gene_degree_stats = args.gene_degree_stats,
            plot_net = args.plot_net,
            gephi_tsv = args.gephi_tsv,
            cytoscape_gml = args.cytoscape_gml,
            map_path = args.map_path)
    elif args.command == "plot":
        plot(
            pickle_path=args.pickle_path,
            with_labels=args.with_labels,
            pdf_path=args.pdf_path,
            gephi_path=args.gephi_path,
            cytoscape_path=args.cytoscape_path,
            max_nodes=args.max_nodes,
            max_edges=args.max_edges,
            symbol=args.symbol,
            map_path=args.map_path,
            species=args.species
        )
    elif args.command == "stats":
       stats(
            dat_file=args.dat_file,
            rewire_net=args.rewire_net,
            out_file=args.out_file_prefix,
            ppif=args.ppif,
            ddif=args.ddif,
            entrezpfamf=args.entrezpfamf,
            map_path=args.map_path,
            species=args.species
        )

if __name__ == '__main__':
    main()
