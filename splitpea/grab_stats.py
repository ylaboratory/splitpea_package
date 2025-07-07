# splitpea/grab_stats.py

import networkx as nx
import pandas as pd
import numpy as np 

def rewired_genes(rewired, background, mapping_path):

    mapping = pd.read_csv(mapping_path, sep="\t", dtype=str)
    entrez_map  = mapping.set_index(mapping['entrez'].astype(str))['symbol'].to_dict()
    ensembl_map = mapping.set_index('ensembl')['symbol'].to_dict()
    uniprot_map = mapping.set_index('uniprot')['symbol'].to_dict()

    data = []
    G = rewired.to_undirected()

    for node in G.nodes():
        node_str = str(node)
        symbol = (
            entrez_map.get(node_str) or
            ensembl_map.get(node_str) or
            uniprot_map.get(node_str) or
            np.nan
        )

        deg    = G.degree(node)
        bg_deg = background.degree(node) if node in background else float('nan')
        norm   = deg / bg_deg if bg_deg > 0 else float('nan')

        gain_count  = 0
        loss_count  = 0
        chaos_count = 0

        for u, v, attr in G.edges(node, data=True):
            if attr.get('chaos', False):
                chaos_count += 1
            else:
                w = attr.get('weight', 0)
                if w > 0:
                    gain_count += 1
                else:
                    loss_count += 1

        if G.has_edge(node, node):
            attr = G[node][node]
            if attr.get('chaos', False):
                chaos_count += 1
            else:
                w = attr.get('weight', 0)
                if w > 0:
                    gain_count += 1
                else:
                    loss_count += 1

        data.append([
            node,
            symbol,
            deg,
            norm,
            gain_count,
            loss_count,
            chaos_count
        ])

    df = pd.DataFrame(
        data,
        columns=[
            "node",
            "symbol",
            "degree",
            "normalized_degree",
            "gain_count",
            "loss_count",
            "chaos_count",
        ]
    )
    return df.sort_values("degree", ascending=False).reset_index(drop=True)

def rewired_edges_stat(dat):
    df = pd.read_csv(dat, sep="\t")
    gain = ((df["weight"] > 0) & (df["chaos"] == False)).sum()
    loss = ((df["weight"] < 0) & (df["chaos"] == False)).sum()
    chaos = (df["chaos"] == True).sum()
    
    return {"gain": gain, "loss": loss, "chaos": chaos}

