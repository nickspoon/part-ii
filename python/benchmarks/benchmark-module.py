import networkx as nx
import numpy as np
import time, sys
import graph, similarity, algebra
from util import *

field = GF(17)
data = []
try:
    for n in range1(2,10,2):
        for d in range1(n-1):
            for i in range(4):
                print (n, d, i)
                G = nx.random_regular_graph(d, n)
                Gp = graph.shuffle_graph(G)

                [B1, B2] = graph.stable_colouring([G, Gp])
                C1 = [numpy_to_nzmath(B1[i], field) for i in range(len(B1)) if np.sum(B1[i]) != 0 ]
                C2 = [numpy_to_nzmath(B2[i], field) for i in range(len(B2)) if np.sum(B1[i]) != 0 ]

                (Lbasis, Vbasis) = similarity.compute_basis3(C1, C2)
                L = algebra.Algebra.fromBasis(field, Lbasis)
                V = algebra.Module.fromBasis(field, Lbasis, Vbasis, L)

                t = time.time()
                v = V.findGenerator()
                t = time.time() - t
                data.append((n, d, L.dim, t))
except (Exception, KeyboardInterrupt):
    if not data:
        sys.exit(0)
    pass

f = open("generator-17.dat", "a")
for entry in data:
    print >>f, "%d %d %d %f" % entry
f.close()

