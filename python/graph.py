import sys, os
import networkx as nx
import numpy as np
from subprocess import Popen, PIPE
import tempfile
from util import *
from similarity import similarity

STABCOL_PATH="../stabprogs/stabcol"

# Load the nth graph from a graph6 file
def load_graph(fn, n):
    G = None
    f = open(fn, "r")
    for i, line in enumerate(f):
        if i == n:
            G = nx.parse_graph6(line.strip())
            break
    if G is None:
        raise EOFError("Graph %d not found: EOF reached" % n)
    f.close()
    return G

# Write a graph in a form STABCOL can process    
def write_stabcol_input(G):
    A = nx.to_numpy_matrix(G, None, int)
    # Set the colour of elements in the diagonal to 2
    # This is required for STABCOL to process the matrix
    np.fill_diagonal(A, 2)
    
    (fd, path) = tempfile.mkstemp()
    f = os.fdopen(fd, "w")
    f.write("3\n") # number of colours
    f.write(str(nx.number_of_nodes(G)) + "\n") # number of vertices
    np.savetxt(f, A, '%d')
    f.close()
    return path

def run_stabcol(fn):
    stabcol = Popen([STABCOL_PATH, fn], stdout=PIPE)
    # Parse stabcol output
    adjmat = False
    for line in stabcol.stdout:
        if adjmat:
            A = np.loadtxt(stabcol.stdout, dtype=int, ndmin=2)
        if line.startswith("number of colors:"):
            colors = int(line.strip().rpartition(" ")[2])
        if line.startswith("number of cells:"):
            cells = int(line.strip().rpartition(" ")[2])
        if line.startswith("adjacency matrix"):
            adjmat = True
    return (colors, cells, A)

def stabcol(G):
    path = write_stabcol_input(G)
    return run_stabcol(path)
    
# Given a matrix M representing coloured adjacency, generate {0,1}-matrices M_i
# for each colour.
def colour_explode(M, n):
    L = []
    for c in range(n):
        C = (M == c).astype(int)
        L.append(C)
    return L

# Generate a stable colouring for a list of graphs.
# Input: A list of NetworkX graphs
# Output: A list of lists of {0,1} colour matrices
def stable_colouring(graphs):
    union = reduce(nx.disjoint_union, graphs)
    (colours, cells, adj) = stabcol(union)
    # N contains the indices of the starts and ends of diagonal blocks
    N = np.cumsum([0] + [ nx.number_of_nodes(G) for G in graphs ])
    assert (N[-1], N[-1]) == adj.shape
    # Separate the diagonal blocks into a list
    As = [ adj[N[i]:N[i+1], N[i]:N[i+1]] for i in range(len(graphs)) ]
    # Generate {0,1} matrices
    Bs = [ colour_explode(A, colours) for A in As ]
    return Bs

# Test if two graphs G1 and G2 are isomorphic using 2-dim WL
def isomorphic(G1, G2):
    if G1.number_of_nodes() != G2.number_of_nodes():
        print "Non-isomorphic: different number of nodes."
        return False
    [B1, B2] = stable_colouring([G1, G2])
    print "%d colours found." % len(B1)
    L1 = [np.sum(M) for M in B1]
    L2 = [np.sum(M) for M in B2]
    if L1 != L2:
        print "Non-isomorphic: different stable colouring."
        return False
    #field = next_prime_field(len(B1))
    field = GF(5)
    p = field.getCharacteristic()
    print "Testing similarity over GF(%d)" % p
    # Convert to numpy matrices, ignoring colours with zero entries
    C1 = [numpy_to_nzmath(B1[i], field) for i in range(len(B1)) if L1[i] != 0 ]
    C2 = [numpy_to_nzmath(B2[i], field) for i in range(len(B2)) if L2[i] != 0 ]
    assert len(C1) == len(C2)
    if similarity(C1, C2) is None:
        print "Non-isomorphic: no simultaneous similarity in GF(%d)." % p
        return False
    print "Isomorphic or too difficult to tell."
    return True

if __name__ == '__main__':
    Gs = []
    for entry in sys.argv[1:]:
        (path, num) = entry.split(":")
        G = load_graph(path, int(num))
        Gs.append(G)
    print isomorphic(Gs[0], Gs[1])
