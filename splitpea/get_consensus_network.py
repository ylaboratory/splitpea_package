import networkx as nx
import pickle
import time
import glob
import sys, os
from pathlib import Path


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
