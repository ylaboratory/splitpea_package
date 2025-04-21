# splitpea/grab_stats.py

import matplotlib.pyplot as plt
import networkx as nx
from adjustText import adjust_text

import matplotlib.pyplot as plt
import networkx as nx
from adjustText import adjust_text

def plot_splice_network(
    G,
    layout=None,
    node_size=30,
    with_labels=False,
    pdf_path=None,
    plot_matplot=True,
    gephi_path=None
):
    """
    Plots G with:
      - blue edges for weight>0 & chaos==False
      - red edges   for weight<0 & chaos==False
      - yellow edges for chaos==True

    Args:
      G            : networkx.Graph
      layout       : dict of positions {node:(x,y)} or None to compute spring_layout
      node_size    : size of nodes
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
        nx.draw_networkx_nodes(G, pos, node_size=node_size)
        nx.draw_networkx_edges(G, pos, edgelist=gain_edges,  edge_color="blue",   label="gain")
        nx.draw_networkx_edges(G, pos, edgelist=loss_edges,  edge_color="red",    label="loss")
        nx.draw_networkx_edges(G, pos, edgelist=chaos_edges, edge_color="yellow", label="chaos")

        if with_labels:
            texts = []
            for n, (x, y) in pos.items():
                texts.append(plt.text(x, y, str(n), fontsize=8))
            adjust_text(texts, arrowprops=dict(arrowstyle="->", color='grey'))

        plt.axis("off")
        plt.tight_layout()

        if pdf_path:
            plt.savefig(pdf_path, format="pdf", bbox_inches="tight")

        if plot_matplot:
            plt.show()
        else:
            plt.close()

    if gephi_path:
        nx.write_gexf(G, gephi_path)
