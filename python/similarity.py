from nzmath import matrix, vector, prime
from linalg import *
from util import *
from itertools import product
from multiprocessing import Pool
from math import ceil
import algebra
import pickle

PROCESSES=8

def pack_matrix(M):
    return (nzmath_to_numpy(M), M.coeff_ring.getCharacteristic())

def unpack_matrix((X, p)):
    return numpy_to_nzmath(X, GF(p))

def weaksim((Ap, Bp)):
    try:
        A = unpack_matrix(Ap)
        if Bp is None: B = A
        else: B = unpack_matrix(Bp)
        lineq = MatrixLinearEquations(A.coeff_ring, A.row, A.column)
        lineq.weakSimilarity(A, B)
        return [ pack_matrix(M) for M in lineq.kernel() ]
    except KeyboardInterrupt:
        pass

def chunks(l, n):
    return (l[i:i+n] for i in range(0, len(l), int(ceil(len(l)/float(n)))))

def parallel_intersect(Lp):
    L = [ map(unpack_matrix, X) for X in Lp ]
    R = reduce(intersect_solutions, L)
    return map(pack_matrix, R)

#def weaksim(A, B=None):
#    try:
#        return A
#    except KeyboardInterrupt:
#        pass

# Given two lists of matrices {A_1, ... A_n}, {B_1, ..., B_n}, find X
# such that for all i = 1..n, X * A_i * X^-1 = B_i 
def similarity(As, Bs, field):
    pool = Pool(processes=PROCESSES)

    dim = As[0].row
    for A in As: assert A.row == A.column == dim
    for B in Bs: assert B.row == B.column == dim
    assert len(As) == len(Bs)
    print "Simultaneous similarity over %d %dx%d matrices" % (len(As), dim, dim)
    Lspaces = pool.map_async(weaksim, [ (pack_matrix(B), None) for B in Bs ])
    Vspaces = pool.map_async(weaksim, 
            [ (pack_matrix(A), pack_matrix(B)) for (A, B) in zip(As, Bs) ])
    spaces = [ map(unpack_matrix, S) for S in Lspaces.get() ]
    # First compute the basis of the algebra L = { X | X*B_i = B_i*X }
    print "Computing the basis of algebra L = { X | X*B_i = B_i*X }"
    #spaces = []
    #for B in Bs:
    #    lineq = MatrixLinearEquations(field, dim, dim)
    #    lineq.weakSimilarity(B)
    #    spaces.append(lineq.kernel())
    spaces_div = chunks([map(pack_matrix, L) for L in spaces], PROCESSES)
    Lspaces = pool.map_async(parallel_intersect, spaces_div)
    
    #lineq = MatrixLinearEquations(field, dim, dim)
    #for B in Bs:
    #    lineq.weakSimilarity(B)
    #Lbasis = lineq.kernel()
    
    # Next compute the basis of the L-module V = { X | X*A_i = B_i*X }
    spaces = [ map(unpack_matrix, S) for S in Vspaces.get() ]
    spaces_div = chunks([map(pack_matrix, L) for L in spaces], PROCESSES)
    Vspaces = pool.map_async(parallel_intersect, spaces_div)
    print "Computing the basis of L-module V = { X | X*A_i = B_i*X }"
#    spaces = []
#    for (A, B) in zip(As, Bs):
#        lineq = MatrixLinearEquations(field, dim, dim)
#        lineq.weakSimilarity(A, B)
#        spaces.append(lineq.kernel())
    result = [ map(unpack_matrix, L) for L in Lspaces.get() ]
    Lbasis = reduce(intersect_solutions, result)
    result = [ map(unpack_matrix, L) for L in Vspaces.get() ]
    Vbasis = reduce(intersect_solutions, result)
    pool.close()
    
    #lineq = MatrixLinearEquations(field, dim, dim)
    #for (A, B) in zip(As, Bs):
    #    lineq.weakSimilarity(A, B)
    #Vbasis = lineq.kernel()
    # If V = {}, then there is no solution.
    if not Vbasis:
        print "No X found such that X*A_i = B_i*X"
        return None
    
    assert all([X * B == B * X for B in Bs for X in Lbasis])
    assert all([X * A == B * X for (A, B) in zip(As, Bs) for X in Vbasis])
    
    print "Constructing %d-dimensional algebra L" % len(Lbasis)
    L = algebra.Algebra.fromBasis(field, Lbasis)
    print "Constructing %d-dimensional module V" % len(Vbasis)
    V = algebra.Module.fromBasis(field, Lbasis, Vbasis, L)
    #print_matrices(V.stconsts)
    print "Finding a generator in V"
    v = V.findGenerator()
    
    if v is not None:
        Z = vector_to_matrix(v, Vbasis)
        # Check if Z is invertible, i.e. if v is a unit. If v is not a unit
        # then there are no units in V.
        try:
            Zprime = Z.inverse()
        except matrix.NoInverse:
            return None
        assert all(Z * A * Zprime == B for (A,B) in zip(As, Bs))
        return Z
    return None

if __name__ == "__main__":
    print "Random matrix test"
    R = random_matrix(3, GF(3))
    try:
        print_matrices([R, R.inverse()])
    except matrix.NoInverse:
        print "Not invertible"

    print "Simultaneous similarity test"
    dim = 12
    n = 5
    field = GF(prime.randPrime(2))
    invertible = True
    X = random_matrix(dim, field)
    try:
        Xprime = X.inverse()
    except matrix.NoInverse:
        invertible = False
    Ms = [ random_matrix(dim, field) for i in range(n) ]
    As = [ X * M for M in Ms ]
    Bs = [ M * X for M in Ms ]
    print X
    print "invertible =", invertible
    print
    print similarity(As, Bs, field)
