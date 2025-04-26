# splitpea/core.py

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
import shutil


def tb_query(tb_file, chrom, start, end):
    """
    Call tabix and yield an array of strings for each line returned.
    (Adapted from https://github.com/slowkow/pytabix)
    """
    if shutil.which("tabix") is None:
        raise RuntimeError("The 'tabix' binary is required but was not found. Please install it via 'sudo apt-get install tabix' or 'brew install htslib'.")
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
