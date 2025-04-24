# splitpea/grab_stats.py

import matplotlib.pyplot as plt
import networkx as nx
from adjustText import adjust_text

import matplotlib.pyplot as plt
import networkx as nx
from adjustText import adjust_text


def simplify_graph(G):
    H = nx.Graph()
    H.add_nodes_from(G.nodes()) 
    
    for u, v, d in G.edges(data=True):
        weight = d.get('weight', 0)
        chaos  = bool(d.get('chaos', False))
        
        if chaos:
            etype = 'chaos'
        elif weight > 0:
            etype = 'gain'
        elif weight < 0:
            etype = 'loss'
        else:
            etype = 'zero'  
        
        H.add_edge(u, v, weight=weight, chaos=chaos, type=etype)
    
    return H


def plot_rewired_network(
    G,
    layout=None,
    node_size=None,
    edge_width=0.25,
    edge_alpha=0.5,            
    with_labels=False,
    pdf_path=None,
    plot_matplot=True,
    gephi_path=None
):
    """
    Plots G using matplotlib and gives option to save a Gephi file for further 
    modifications. 

    Args:
      G            : networkx.Graph
      layout       : dict of positions {node:(x,y)} or None to compute spring_layout
      node_size    : size of nodes
      edge_width   : width of edges 
      edge_alpha   : opacity (0.0 transparent → 1.0 opaque)
      with_labels  : whether to draw node labels 
      pdf_path     : if str, save the plot as a PDF to this path
      plot_matplot : whether to plot using Matplotlib
      gephi_path   : if str, export G as GEXF to this path (for Gephi)
    """

    if not plot_matplot and pdf_path is None and gephi_path is None:
        raise ValueError("Must either plot with Matplotlib (plot_matplot=True), save as PDF (pdf_path), or export to Gephi (gephi_path)")

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

    if plot_matplot or pdf_path:
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
        else:
            plt.close()

    if gephi_path:
        simpleG = simplify_graph(G)
        nx.write_gexf(simpleG, gephi_path)

        