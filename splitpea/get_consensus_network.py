import networkx as nx
import pickle
import time
import glob
import sys, os
from pathlib import Path
from typing import Iterable, Optional
import pandas as pd
import matplotlib.pyplot as plt


def get_consensus_network(net_dir):
    """
    net_dir: an input directory containing sample-level networks
    """

    print("[" + time.strftime("%H:%M:%S", time.localtime()) + "] Loading networks...")

    c_consensus_neg = nx.Graph(name="consensus_neg")
    c_consensus_neg.graph["num_graphs"] = 0

    c_consensus_pos = nx.Graph(name="consensus_pos")
    c_consensus_pos.graph["num_graphs"] = 0

    for f in glob.glob(net_dir + "/*.pickle"):
        print(
            "[" + time.strftime("%H:%M:%S", time.localtime()) + "] Adding " + f + "...."
        )
        g = pickle.load(open(f, "rb"))

        g_neg = g.copy()
        g_pos = g.copy()

        for gi, gj in g.edges():
            if g[gi][gj]["chaos"]:
                g_neg.remove_edge(gi, gj)
                g_pos.remove_edge(gi, gj)
                continue

            if g[gi][gj]["weight"] <= 0:
                g_pos.remove_edge(gi, gj)
                g_neg[gi][gj]["num_neg"] = 1
            else:
                g_neg.remove_edge(gi, gj)
                g_pos[gi][gj]["num_pos"] = 1

        g_neg.remove_nodes_from([n for (n, deg) in g_neg.degree() if deg == 0])
        g_pos.remove_nodes_from([n for (n, deg) in g_pos.degree() if deg == 0])

        comb_neg = nx.compose(g_neg, c_consensus_neg)
        neg_edge_weight = {
            e: c_consensus_neg.edges[e]["weight"] + g_neg.edges[e]["weight"]
            for e in c_consensus_neg.edges & g_neg.edges
        }
        neg_edge_num = {
            e: c_consensus_neg.edges[e]["num_neg"] + g_neg.edges[e]["num_neg"]
            for e in c_consensus_neg.edges & g_neg.edges
        }
        nx.set_edge_attributes(comb_neg, neg_edge_weight, "weight")
        nx.set_edge_attributes(comb_neg, neg_edge_num, "num_neg")
        c_consensus_neg = comb_neg.copy()
        c_consensus_neg.graph["num_graphs"] += 1

        comb_pos = nx.compose(g_pos, c_consensus_pos)
        pos_edge_weight = {
            e: c_consensus_pos.edges[e]["weight"] + g_pos.edges[e]["weight"]
            for e in c_consensus_pos.edges & g_pos.edges
        }
        pos_edge_num = {
            e: c_consensus_pos.edges[e]["num_pos"] + g_pos.edges[e]["num_pos"]
            for e in c_consensus_pos.edges & g_pos.edges
        }
        nx.set_edge_attributes(comb_pos, pos_edge_weight, "weight")
        nx.set_edge_attributes(comb_pos, pos_edge_num, "num_pos")
        c_consensus_pos = comb_pos.copy()
        c_consensus_pos.graph["num_graphs"] += 1

    with open("consensus_network_neg.pickle", "wb") as out:
        pickle.dump(c_consensus_neg, out)

    with open("consensus_network_pos.pickle", "wb") as out:
        pickle.dump(c_consensus_pos, out)

    return c_consensus_neg, c_consensus_pos


def analyze_consensus_threshold(
    neg_path: Optional[str] = None,
    pos_path: Optional[str] = None,
    neg_networkx: Optional[nx.Graph] = None,
    pos_networkx: Optional[nx.Graph] = None,
    thresholds: Iterable[float] = (0.1, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0),
    label: str = "sample",
    # outputs
    pickles_dir: Optional[str] = None,  # default: ./threshold_networks/<label>/
    save_pickles: bool = True,
    write_txt: bool = True,
    txt_path: Optional[
        str
    ] = None,  # default: <pickles_dir>/consensus_threshold.lcc_sizes.txt
    # plotting
    save_pdf_prefix: Optional[
        str
    ] = None,  # saves ..._num_nodes.pdf and ..._prop_nodes.pdf
    title_prefix: Optional[str] = None,  # adds "<title_prefix>-" to plot titles
) -> pd.DataFrame:
    """
    Analyze consensus thresholds for one label (e.g., 'BRCA').

    Provide either/both:
      - neg_path: path to <label>_consensus_neg.pickle (edge attr 'num_neg')
      - pos_path: path to <label>_consensus_pos.pickle (edge attr 'num_pos')

    Saves per-threshold graphs (no joint) into pickles_dir/<direction>/,
    writes a TSV of sizes, and plots BOTH #nodes and proportion-of-nodes
    (relative to threshold=0 baseline).

    Returns a DataFrame with columns:
      [label, direction, threshold, num_nodes, num_edges, prop_nodes]
    """
    if not neg_path and not pos_path and not neg_networkx and not pos_networkx:
        raise ValueError(
            "Provide at least one of neg_path, pos_path, neg_networkx, or pos_networkx."
        )

    def _largest_cc_or_empty(G: nx.Graph) -> nx.Graph:
        if G is None or G.number_of_edges() == 0 or G.number_of_nodes() == 0:
            return nx.Graph()
        H = G.copy()
        H.remove_nodes_from([n for n, d in H.degree() if d == 0])
        if H.number_of_edges() == 0 or H.number_of_nodes() == 0:
            return nx.Graph()
        return H.subgraph(max(nx.connected_components(H), key=len)).copy()

    def _threshold_graph(
        G: nx.Graph, edge_count_attr: str, num_graphs: int, thr: float
    ) -> nx.Graph:
        if G is None or G.number_of_edges() == 0:
            return nx.Graph()
        cutoff = thr * num_graphs
        keep_edges = [
            (u, v)
            for u, v in G.edges
            if G.edges[(u, v)].get(edge_count_attr, 0) >= cutoff
        ]
        if not keep_edges:
            return nx.Graph()
        H = nx.Graph()
        H.add_nodes_from(G.nodes(data=True))
        for u, v in keep_edges:
            H.add_edge(u, v, **G.edges[(u, v)])
        return _largest_cc_or_empty(H)

    graphs = {}
    edge_attr = {}
    if neg_path or neg_networkx:
        if neg_networkx:
            graphs["negative"] = neg_networkx
        else:
            with open(neg_path, "rb") as f:
                graphs["negative"] = pickle.load(f)
        edge_attr["negative"] = "num_neg"
    if pos_path or pos_networkx:
        if pos_networkx:
            graphs["positive"] = pos_networkx
        else:
            with open(pos_path, "rb") as f:
                graphs["positive"] = pickle.load(f)
        edge_attr["positive"] = "num_pos"

    nums = []
    for d, G in graphs.items():
        ng = G.graph.get("num_graphs")
        if ng is None:
            raise ValueError(f"'num_graphs' missing on {d} graph.")
        nums.append(int(ng))
    if len(nums) > 1 and len(set(nums)) != 1:
        raise ValueError("Expected identical 'num_graphs' across provided graphs.")
    num_graphs = nums[0]

    # where to save pickles
    if pickles_dir is None:
        pickles_dir = str(Path.cwd() / "threshold_networks" / label)
    os.makedirs(pickles_dir, exist_ok=True)

    # TSV path default
    if write_txt and txt_path is None:
        txt_path = os.path.join(pickles_dir, "consensus_threshold.lcc_sizes.txt")

    rows = []
    baseline_nodes = {}

    for d, G in graphs.items():
        n0 = G.number_of_nodes()
        baseline_nodes[d] = max(n0, 1)
        rows.append(
            {
                "label": label,
                "direction": d,
                "threshold": 0.0,
                "num_nodes": n0,
                "num_edges": G.number_of_edges(),
                "prop_nodes": n0 / baseline_nodes[d],
            }
        )

    for t in thresholds:
        for d, G in graphs.items():
            H = _threshold_graph(G, edge_attr[d], num_graphs, t)

            if save_pickles:
                d_dir = os.path.join(pickles_dir, d)
                os.makedirs(d_dir, exist_ok=True)
                out_pickle = os.path.join(
                    d_dir,
                    f"{label}_consensus_{'neg' if d=='negative' else 'pos'}.thres_{t}.pickle",
                )
                with open(out_pickle, "wb") as f:
                    pickle.dump(H, f)

            n = H.number_of_nodes()
            rows.append(
                {
                    "label": label,
                    "direction": d,
                    "threshold": float(t),
                    "num_nodes": n,
                    "num_edges": H.number_of_edges(),
                    "prop_nodes": n / baseline_nodes[d],
                }
            )

    df = pd.DataFrame(
        rows,
        columns=[
            "label",
            "direction",
            "threshold",
            "num_nodes",
            "num_edges",
            "prop_nodes",
        ],
    )

    if write_txt:
        header = not os.path.exists(txt_path)
        df.to_csv(
            txt_path,
            sep="\t",
            index=False,
            mode="a" if not header else "w",
            header=header,
        )

    def _plot(df_plot: pd.DataFrame, y: str, title_suffix: str, pdf_suffix: str):
        plt.figure(figsize=(7.5, 5))
        for d in sorted(df_plot["direction"].unique()):
            dd = df_plot[df_plot["direction"] == d].sort_values("threshold")
            plt.plot(dd["threshold"], dd[y], marker="o", label=d)
        plt.xlabel("consensus threshold")
        plt.ylabel("# nodes" if y == "num_nodes" else "proportion of nodes")
        plt.title((title_prefix + " - " if title_prefix else "") + f"{title_suffix}")
        plt.legend(title="direction")
        plt.tight_layout()
        if save_pdf_prefix:
            out_pdf = f"{save_pdf_prefix}_{pdf_suffix}.pdf"
            plt.savefig(out_pdf, bbox_inches="tight")
        plt.show()

    _plot(df, "num_nodes", "nodes vs threshold", "num_nodes")
    _plot(df, "prop_nodes", "proportion of nodes vs threshold", "prop_nodes")

    return df
