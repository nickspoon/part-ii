from nzmath import matrix, vector, prime
from linalg import *
from util import *
from itertools import product
import algebra
import pickle
import pool

def weaksim_packed_list(Aps, Bps):
    As = map(unpack_matrix, Aps)
    if Bps is None: Bs = As
    else: Bs = map(unpack_matrix, Bps)
    dim = As[0].row
    field = As[0].coeff_ring
    lineq = MatrixLinearEquations(field, dim, dim)
    for (A, B) in zip(As, Bs):
        lineq.weakSimilarity(A, B)
    result = lineq.kernel()
    return map(pack_matrix, result)

def compute_basis1(As, Bs):
    dim = As[0].row
    field = As[0].coeff_ring
    # First compute the basis of the algebra L = { X | X*B_i = B_i*X }
    lineq = MatrixLinearEquations(field, dim, dim)
    for B in Bs:
        lineq.weakSimilarity(B)
    Lbasis = lineq.kernel()
    # Next compute the basis of the L-module V = { X | X*A_i = B_i*X }
    lineq = MatrixLinearEquations(field, dim, dim)
    for (A, B) in zip(As, Bs):
        lineq.weakSimilarity(A, B)
    Vbasis = lineq.kernel()
    return (Lbasis, Vbasis)

def compute_basis2(As, Bs):
    Lspaces = parallel_weaksim(Bs)
    Vspaces = parallel_weaksim(As, Bs)
    Lbasis = ParallelIntersection(Lspaces.get(), packed=True, nch=pool.PROCESSES/2)
    Vbasis = ParallelIntersection(Vspaces.get(), packed=True, nch=pool.PROCESSES/2)
    return (Lbasis.get(), Vbasis.get())

def compute_basis3(As, Bs):
    Lspaces = pool.pool().apply_async(weaksim_packed_list,
                    [map(pack_matrix, Bs), None])
    Vspaces = pool.pool().apply_async(weaksim_packed_list,
        [ map(pack_matrix, As), map(pack_matrix, Bs) ])
    Lbasis = map(unpack_matrix, Lspaces.get())
    Vbasis = map(unpack_matrix, Vspaces.get())
    return (Lbasis, Vbasis)

# Given two lists of matrices {A_1, ... A_n}, {B_1, ..., B_n}, find X
# such that for all i = 1..n, X * A_i * X^-1 = B_i 
def similarity(As, Bs):
    pool.start_pool()

    dim = As[0].row
    field = As[0].coeff_ring
    for A in As: assert A.row == A.column == dim
    for B in Bs: assert B.row == B.column == dim
    assert len(As) == len(Bs)
    print "Simultaneous similarity over %d %dx%d matrices" % (len(As), dim, dim)
    (Lbasis, Vbasis) = compute_basis2(As, Bs)
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
    dim = 4
    n = 2
    field = GF(2)
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
    print similarity(As, Bs)
