# splitpea/grab_stats.py

import matplotlib.pyplot as plt
import networkx as nx
from adjustText import adjust_text
import numpy as np

def return_gephi(G, outfile_path):

    with open(outfile_path, 'w') as f_out:
        header = ["source", "target", "weight", "direction"]
        f_out.write('\t'.join(header) + '\n')
        
        for gi, gj in nx.edges(G):
            weight = G[gi][gj]['weight']
            direction = "positive" if weight >= 0 else "negative"
            
            if G[gi][gj].get('chaos', False):
                direction = "chaos-" + direction

            ppi_ddi_out = [gi, gj, str(abs(weight)), direction]
            f_out.write('\t'.join(ppi_ddi_out) + '\n')

def return_cytoscape(G, outfile_path):

    for n, data in G.nodes(data=True):
        data.clear()

    for u, v, data in G.edges(data=True):

        if 'ppi_ddis' in data:
            data.pop('ppi_ddis')

        for key in list(data.keys()):
            if '_' in key:
                new_key = key.replace('_', '')
                data[new_key] = data.pop(key)

        for key in list(data.keys()):
            value = data[key]
            if isinstance(value, (np.integer, np.int32, np.int64)):
                data[key] = int(value)
            elif isinstance(value, (np.float32, np.float64)):
                data[key] = float(value)
            elif isinstance(value, np.bool_):
                data[key] = bool(value)

    nx.write_gml(G, outfile_path)


def plot_rewired_network(
    G,
    layout=None,
    node_size=None,
    edge_width=0.25,
    edge_alpha=0.5,            
    with_labels=False,
    pdf_path=None,
    plot_matplot=True,
    gephi_path=None,
    cytoscape_path=None,
    max_nodes=2000,
    max_edges=10000
):
    """
    Plots G using matplotlib and gives option to save a Gephi file for further 
    modifications. Warns if the graph exceeds max_nodes or max_edges and disables
    matplotlib plotting in that case.

    Args:
      G             : networkx.Graph
      layout        : dict of positions {node:(x,y)} or None to compute spring_layout
      node_size     : size of nodes
      edge_width    : width of edges 
      edge_alpha    : opacity (0.0 transparent → 1.0 opaque)
      with_labels   : whether to draw node labels 
      pdf_path      : if str, save the plot as a PDF to this path
      plot_matplot  : whether to plot using Matplotlib
      gephi_path    : if str, export G as GEXF to this path (for Gephi)
      cytoscape_path: if str, export G as GML to this path (for Cytoscape)
      max_nodes     : maximum node count before disabling matplotlib plotting
      max_edges     : maximum edge count before disabling matplotlib plotting
    """

    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    over_nodes = n_nodes > max_nodes
    over_edges = n_edges > max_edges

    if over_nodes or over_edges:
        parts = []
        if over_nodes:
            parts.append(f"{n_nodes} nodes (limit={max_nodes})")
        if over_edges:
            parts.append(f"{n_edges} edges (limit={max_edges})")
        detail = " and ".join(parts)
        print(
            f"WARNING: graph is large ({detail}). "
            "Plotting may be slow or use too much RAM. "
            "Consider exporting to Gephi/Cytoscape, "
            "or increase max_nodes/max_edges."
        )
        plot_matplot = False

    if not plot_matplot and pdf_path is None and gephi_path is None and cytoscape_path is None:
        raise ValueError(
            "Must either plot with Matplotlib (plot_matplot=True), save as PDF (pdf_path), "
            "or export to Gephi (gephi_path) or Cytoscape (cytoscape_path)"
        )

    gain_edges  = [(u, v) for u, v, d in G.edges(data=True)
                   if d.get("weight", 0) > 0 and not d.get("chaos", False)]
    loss_edges  = [(u, v) for u, v, d in G.edges(data=True)
                   if d.get("weight", 0) < 0 and not d.get("chaos", False)]
    chaos_edges = [(u, v) for u, v, d in G.edges(data=True)
                   if d.get("chaos", False)]

    if layout is None:
        pos = nx.spring_layout(G)
    else:
        pos = layout

    if plot_matplot and pdf_path:
        plt.figure()

        if node_size is None:
            min_size, max_size = 5, 300
            degree_dict = dict(G.degree())
            deg_vals = degree_dict.values()
            scale = lambda d: min_size + (d / max(deg_vals)) * (max_size - min_size)
            node_size = [scale(degree_dict[n]) for n in G.nodes()]

        nx.draw_networkx_nodes(G, pos, node_size=node_size)

        nx.draw_networkx_edges(
            G, pos, edgelist=gain_edges,
            edge_color="blue", width=edge_width,
            alpha=edge_alpha, label="gain"
        )
        nx.draw_networkx_edges(
            G, pos, edgelist=loss_edges,
            edge_color="red", width=edge_width,
            alpha=edge_alpha, label="loss"
        )
        nx.draw_networkx_edges(
            G, pos, edgelist=chaos_edges,
            edge_color="yellow", width=edge_width,
            alpha=edge_alpha, label="chaos"
        )

        if with_labels:
            texts = []
            for n, (x, y) in pos.items():
                texts.append(plt.text(x, y, str(n), fontsize=8))
            adjust_text(texts, arrowprops=dict(arrowstyle="->", color="grey"))

        plt.axis("off")
        plt.tight_layout()

        if pdf_path:
            plt.savefig(pdf_path, format="pdf", bbox_inches="tight")

        if plot_matplot:
            plt.show()
            plt.close()
        else:
            plt.close()

    if gephi_path:
        _ = return_gephi(G, gephi_path)
    
    if cytoscape_path:
        _ = return_cytoscape(G, cytoscape_path)
