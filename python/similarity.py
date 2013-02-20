from nzmath import matrix, vector, prime
from linalg import *
from util import *
from itertools import product
import algebra

# Given two lists of matrices {A_1, ... A_n}, {B_1, ..., B_n}, find X
# such that for all i = 1..n, X * A_i * X^-1 = B_i 
def similarity(As, Bs, field):
    dim = As[0].row
    for A in As: assert A.row == A.column == dim
    for B in Bs: assert B.row == B.column == dim
    assert len(As) == len(Bs)
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
    # If V = {}, then there is no solution.
    if not Vbasis: return None
    
    assert all([X * B == B * X for B in Bs for X in Lbasis])
    assert all([X * A == B * X for (A, B) in zip(As, Bs) for X in Vbasis])
    
    L = algebra.Algebra.fromBasis(field, Lbasis)
    V = algebra.Module.fromBasis(field, Lbasis, Vbasis, L)
    #print_matrices(V.stconsts)
    v = V.findGenerator()
    
    if v is not None:
        Z = vector_to_matrix(v, Vbasis)
        # Check if Z is invertible, i.e. if v is a unit. If v is not a unit
        # then there are no units in V.
        try:
            Zprime = Z.inverse()
        except matrix.NoInverse:
            return None
        assert all(Z * A * Z.inverse() == B for (A,B) in zip(As, Bs))
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
    dim = 7
    n = 4
    field = GF(prime.randPrime(3))
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
