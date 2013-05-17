import graph
import itertools, sys, random
import networkx as nx
import matplotlib.pyplot as plt

def test_graphs(fn):
    graphs = graph.load_all_graphs(fn)
    for (i,j) in itertools.product(range(len(graphs)), repeat=2):
        iso = graph.isomorphic(graphs[i], graphs[j])
        assert iso == (i == j), "Incorrect result for " + str(i) + "," + str(j)

def test_iso(n):
    G = nx.gnp_random_graph(n, 0.4)
    Gp = graph.shuffle_graph(G)
    assert graph.isomorphic(G, Gp)
    #U = nx.disjoint_union(G, Gp)
    #pos=nx.graphviz_layout(U,prog="neato")
    #nx.draw(G, nx.circular_layout(G), node_size=1000, font_size=18)
    #plt.savefig("isograph1.eps")
    #plt.figure()
    #nx.draw(G, nx.spring_layout(G), node_size=1000, font_size=18)
    #plt.savefig("isograph2.eps")

if __name__ == "__main__":
    test_graphs(sys.argv[1])
    #test_iso(int(sys.argv[1]))
