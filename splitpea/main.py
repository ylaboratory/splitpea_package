import argparse
import pickle
from collections import defaultdict
import networkx as nx
import os
from numpy import average
from functools import partial
from os.path import basename
from subprocess import Popen, PIPE
from pathlib import Path
import logging
import csv
import tempfile

from .exons import Exons
from .parser import parse_suppa2
from .parser import parse_rmats



logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s',
                    datefmt="%m-%d-%Y %I:%M:%S %p",
                    level=logging.ERROR)
logger = logging.getLogger('diff_exon')


def tb_query(tb_file, chrom, start, end):
    """
    Call tabix and yield an array of strings for each line returned.
    (Adapted from https://github.com/slowkow/pytabix)
    """
    query = '{}:{}-{}'.format(chrom, start, end)
    process = Popen(['tabix', '-f', tb_file, query], stdout=PIPE, text=True)
    for line in process.stdout:
        yield line.strip().split()


def calculate_edges(g, dex, tb, all_ddis, all_ppis, entrez_pfams):
    """
    For each differentially expressed exon, build a network based on reference PPI and DDI.
    """
    for ex in dex.dexon2scores:
        ex_dpsi = dex.dexon2scores[ex][0]

        for r in tb_query(tb, ex[0], ex[1], ex[2]):
            if r[8] != ex[3]:  # Ensure same strand
                continue

            gi = r[0]
            pfi = r[2]
            pfi_start = r[3]
            pfi_end = r[4]

            if gi not in all_ppis or gi not in entrez_pfams:
                continue

            if pfi not in all_ddis or pfi not in entrez_pfams:
                continue

            for gj in all_ppis[gi]:
                if gj not in entrez_pfams:
                    continue

                ddi_overlap = set(entrez_pfams[gj]).intersection(set(all_ddis[pfi]))
                if len(ddi_overlap) > 0:
                    if gi not in g:
                        g.add_node(gi)

                    if pfi not in g.nodes[gi]:
                        g.nodes[gi][pfi] = defaultdict(dict)

                    if (pfi_start, pfi_end) not in g.nodes[gi][pfi]:
                        g.nodes[gi][pfi][(pfi_start, pfi_end)]['dexons'] = set()
                        g.nodes[gi][pfi][(pfi_start, pfi_end)]['min_dpsi'] = 1

                    g.nodes[gi][pfi][(pfi_start, pfi_end)]['dexons'].add((ex, ex_dpsi))
                    if ex_dpsi < g.nodes[gi][pfi][(pfi_start, pfi_end)]['min_dpsi']:
                        g.nodes[gi][pfi][(pfi_start, pfi_end)]['min_dpsi'] = ex_dpsi

                    if gj not in g[gi]:
                        g.add_edge(gi, gj,
                                   ppi_weight=all_ppis[gi][gj]['weight'],
                                   ppi_ddis=defaultdict(partial(defaultdict, set)))

                    g.edges[gi, gj]['ppi_ddis'][gi][pfi] |= ddi_overlap


def deduce_final_edge_weights(g):
    """
    Collapse multiple domain edges into final edge weights.
    """
    num_chaos = 0
    for gi, gj in g.edges():
        chaos = False
        pos_dpsi = []
        neg_dpsi = []

        for g_ddi in g[gi][gj]['ppi_ddis']:
            for pf_g in g[gi][gj]['ppi_ddis'][g_ddi]:
                for pf_loc in g.nodes[g_ddi][pf_g]:
                    pf_dpsi = g.nodes[g_ddi][pf_g][pf_loc]['min_dpsi']
                    if pf_dpsi < 0:
                        neg_dpsi.append(pf_dpsi)
                    else:
                        pos_dpsi.append(pf_dpsi)

        if pos_dpsi and neg_dpsi:
            chaos = True
            num_chaos += 1

        g[gi][gj]['weight'] = average(pos_dpsi + neg_dpsi)
        g[gi][gj]['chaos'] = chaos
        g[gi][gj]['num_dpsi_pfs'] = len(pos_dpsi) + len(neg_dpsi)

    g.graph['num_chaos'] = num_chaos


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
        species: str = 'human'):
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
      
      The following file paths default to files under 'src/reference':
      ppif: Protein-protein interaction file (default: human_ppi_0.5.dat).
      ddif: Domain-domain interaction file (default: ddi_0.5.dat).
      entrezpfamf: Gene-protein domain info file (default: human_entrez_pfam.txt).
      pfamcoordsf: Pfam genome coordinates file (default: human_pfam_genome_coords_sorted.txt.gz).
      tbf: Tabix file (default: human_pfam_genome_coords_sorted.txt.gz).
    """
    if verbose:
        logger.setLevel(logging.INFO)

    if input_format.lower() == "suppa2":
        if not (isinstance(in_file, (list, tuple)) and len(in_file) == 2):
            raise ValueError("For SUPPA2 input format, in_file must be a list or tuple of 2 file paths: [psivec_file, dpsi_file].")
        psivec_file, dpsi_file = in_file
        in_file = parse_suppa2(psivec_file, dpsi_file, species=species)
        skip = 1
    
    if input_format.lower() == "rmats":
        in_file = parse_rmats(in_file)
        skip = 1

    # Set default reference file paths if not provided.
    if ppif is None or ddif is None or entrezpfamf is None or pfamcoordsf is None or tbf is None:
        base_dir = str(Path(__file__).resolve().parent.parent / 'src' / 'reference')
        if ppif is None:
            ppif = os.path.join(base_dir, "human_ppi_0.5.dat")
        if ddif is None:
            ddif = os.path.join(base_dir, "ddi_0.5.dat")
        if entrezpfamf is None:
            entrezpfamf = os.path.join(base_dir, "human_entrez_pfam.txt")
        if pfamcoordsf is None:
            pfamcoordsf = os.path.join(base_dir, "human_pfam_genome_coords_sorted.txt.gz")
        if tbf is None:
            tbf = os.path.join(base_dir, "human_pfam_genome_coords_sorted.txt.gz")

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
    dexons.read_dexons(in_file, skip, dpsi_cut, sigscore_cut, include_nas)

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

    logger.info("Done")
