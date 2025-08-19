# splitpea/get_background_ppi.py

# makes output network w/ scores w/o any diff exon info
# to get what the background network would be

import networkx as nx
from collections import defaultdict
from functools import partial

import time
from pathlib import Path
import pickle
import logging

# # file path constants
# base_dir = str(Path(__file__).resolve().parent.parent)
# ppif = base_dir + "/reference/human_ppi_0.5.dat"
# ddif = base_dir + "/reference/ddi_0.5.dat"
# entrezpfamf = base_dir + "/reference/human_entrez_pfam.txt"

# logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s',
#                     datefmt="%m-%d-%Y %I:%M:%S %p",
#                     level=logging.ERROR)
# # logger = logging.get# logger('diff_exon')


def get_background(ppif, ddif, entrezpfamf):
    # logger.info("Reading in network files for background PPI...")
    all_ddis = nx.read_weighted_edgelist(ddif)
    all_ppis = nx.read_weighted_edgelist(ppif) 
    entrez_pfams = nx.bipartite.read_edgelist(entrezpfamf)
    # logger.info("Full PPI info...")
    # logger.info(all_ppis)

    # logger.info("Calculating background network...")    
    # constructing graph w/ ppi nodes incident on genes w/ a ddi (w/o exon info)
    g = nx.Graph()

    for gi, gj in all_ppis.edges:
        if gi not in entrez_pfams or gj not in entrez_pfams:
            continue
        
        for pfi in entrez_pfams[gi]:
            if pfi not in all_ddis:
                continue

            ddi_overlap = set(entrez_pfams[gj]).intersection(set(all_ddis[pfi]))
            if len(ddi_overlap) > 0:
                if gi not in g:
                    g.add_node(gi)

                if gj not in g[gi]:
                    g.add_edge(gi, gj,
                                ppi_weight = all_ppis[gi][gj]['weight'],
                                ppi_ddis = defaultdict(partial(defaultdict, set))) # cannot use lambda if want to pickle later
                
                g.edges[gi,gj]['ppi_ddis'][gi][pfi] |= ddi_overlap

    # logger.info( "PPI with DDI support (i.e., background)...")
    # logger.info(g)
    # logger.info("Outputting...")

    with open("human_ppi_ddi_bg.pickle", 'wb') as out:
        pickle.dump(g, out)

    return g

# [21:23:07] Reading in network files...
# [21:23:09] Full PPI info...
# Graph with 20286 nodes and 793078 edges
# [21:23:10] Calculating background network...
# [21:23:14] PPI with DDI support (i.e., background)...
# Graph with 10602 nodes and 64739 edges
# [21:23:14] Outputting...
# [21:23:14] Done