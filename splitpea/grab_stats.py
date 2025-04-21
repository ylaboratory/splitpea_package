# splitpea/grab_stats.py

import networkx as nx
import pandas as pd

def rewired_genes(rewired, background):
    data = []
    
    for node in rewired.nodes():
        rewired_degree = rewired.degree(node)
        background_degree = background.degree(node) if node in background else float('nan')  
        normalized_degree = rewired_degree / background_degree if background_degree > 0 else float('nan')
        
        data.append([node, rewired_degree, normalized_degree])
    
    df = pd.DataFrame(data, columns=["node", "degree", "normalized_degree"])
    return df

def rewired_edges_stat(dat):
    df = pd.read_csv(dat, sep="\t")
    gain = ((df["weight"] > 0) & (df["chaos"] == False)).sum()
    loss = ((df["weight"] < 0) & (df["chaos"] == False)).sum()
    chaos = (df["chaos"] == True).sum()
    
    return {"gain": gain, "loss": loss, "chaos": chaos}

