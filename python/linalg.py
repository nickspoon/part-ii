from nzmath import matrix, vector, ring
from util import *
from itertools import product
import pool
from parallel import inject_parallel
import ffpack_interface

# Class representing a system of linear equations.
# The solution is x where self.A * x = self.b
class LinearEquations(object):
    
    # field: Field over which to solve
    # dim: number of variables
    def __init__(self, field, dim, parallel=True):
        self.field = field
        self.dim = dim
        self.A = None
        self.b = None
        self.parallel = parallel
        print "Solving a linear equation system in %d variables" % dim
        
    # X: k*n matrix; y: vector of length k
    def addEquations(self, X, y):
        if X.column != self.dim:
            raise DimensionError("Supplied matrix has the wrong number of columns")
        if X.row != len(y):
            raise DimensionError("Matrix/vector size mismatch")
        if self.A is None:
            self.A = X.copy()
            if self.parallel:
                inject_parallel(self.A)
            self.b = y.copy()
        else:
            self.A.extendRow(X)
            self.b.compo.extend(y)
        
    # x: vector of length n; y: scalar
    def addEquation(self, x, y):
        X = x.toMatrix()
        Y = vector.Vector([y])
        self.addEquations(X, Y)
    
    # For two solution spaces A, B for some matrix X, find coefficients
    # a_i, b_j such that sum_i a_i * A_i = sum_j b_j * B_j
    def intersect(self, S1, S2):
        cvectors = [ flatten_matrix(S) for S in S1 ] + \
                    [ flatten_matrix(-S) for S in S2 ]
        coeffs = matrix.Matrix(len(cvectors[0]), len(cvectors),
                                cvectors, self.field)
        v = flatten_matrix(matrix.zeroMatrix(S1[0].row, S1[0].column))
        self.addEquations(coeffs, v)
    
    def matEqSpace(self, v, S):
        coeffs = matrix.Matrix(len(v), len(S), [ B * v for B in S ])
        self.addEquations(coeffs, v)
    
    # Returns a tuple (solution, kernel)
    def solve(self):
        return self.A.solve(self.b)
        
    # Returns a single solution
    def solution(self):
        #return self.A.inverseImage(self.b)[1]
        return ffpack_interface.solution(self.A, self.b)[1]
    
    def kernel(self):
        #K = self.A.kernel()
        K = ffpack_interface.kernel(self.A)
        if K is None: return None
        return [ K[i] for i in range1(K.column) ]

# As above, but the variables form a matrix.
class MatrixLinearEquations(LinearEquations):
    def __init__(self, field, row, col, parallel=True):
        self.row = row
        self.col = col
        super(MatrixLinearEquations, self).__init__(field, row*col, parallel)
        
    def addCoefficientMatrix(self, A, n):
        self.addEquation(flatten_matrix(A), n)

    def weakSimilarity(self, A, B=None):
        # Find some X s.t. XA = BX, i.e. XA - BX = 0 where X, A are
        # matrices over a field.
        # A, B are square matrices of dimensions m, n
        # Output X has dimension n x m
        if A.column != A.row:
            raise DimensionError("Square matrix required for similarity")
        if B is None:
            B = A
        else:
            if B.column != B.row:
                raise DimensionError("Square matrix required for similarity")
        for i in range1(B.row):
            for j in range1(A.column):
                # Coefficients of X_ij, initially zero
                coeffs = matrix.FieldMatrix(B.row, A.column, self.field)
                # XA: Multiply X-row by A-column
                coeffs.setRow(i, A.getColumn(j))
                # -BX: Multiply X-column by B-row
                coeffs.setColumn(j, -B.getRow(i))
                # XA - BX: handle row/column intersection
                coeffs[i,j] = A[j,j] - B[i,i]
                self.addEquation(flatten_matrix(coeffs), 0)
                
    def matrixEquation(self, a, b):
        # Add a constraint of the form Xa = b, where a and b are vectors
        if len(a) != self.col:
            raise DimensionError("Dimension mismatch: matrix has %d columns \
                        but argument vector has length %d", self.col, len(a))
        if len(b) != self.row:
            raise DimensionError("Dimension mismatch: matrix has %d rows \
                        but result vector has length %d", self.col, len(a))
        for i in range1(self.row):
            coeffs = matrix.FieldSquareMatrix(self.row, self.col, self.field)
            coeffs.setRow(i, a)
            self.addEquation(flatten_matrix(coeffs), b[i])
                
    def solve(self):
        (v, K) = super(MatrixLinearEquations, self).solve()
        return (unflatten_matrix(v, self.row, self.col),
                [unflatten_matrix(w, self.row, self.col)
                    for w in K])
    
    def solution(self):
        v = super(MatrixLinearEquations, self).solution()
        return unflatten_matrix(v[1], self.row, self.col)
    
    def kernel(self):
        K = super(MatrixLinearEquations, self).kernel()
        return [unflatten_matrix(w, self.row, self.col)
                    for w in K]

class ParallelIntersection:
    def __init__(self, spaces, packed=False, nch=pool.PROCESSES):
        if not packed: spaces = [map(pack_matrix, L) for L in spaces]
        self.nch = nch
        spaces_div = list(nchunks(spaces, nch))
        for subl in spaces_div:
            subl.sort(key=len)
        self.result = pool.pool().map_async(intersect_worker, spaces_div)
    
    def get(self):
        spaces = self.result.get()
        if len(spaces) <= 2:
            spc = [ map(unpack_matrix, L) for L in spaces ]
            return reduce(intersect_solutions, spc)
        else:
            p = ParallelIntersection(spaces, packed=True, nch=self.nch/2)
            return p.get()

# Flatten a matrix column-wise into a vector
# e.g. [ [ a b ], [ c d ] ] => [ a c b d ]
def flatten_matrix(A):
    v = vector.Vector([])
    for i in range1(A.column):
        v.compo.extend(A[i])
    return v
    
def unflatten_matrix(v, row, col):
    if len(v) != row * col:
        raise DimensionError("Incorrect matrix dimensions supplied.")
    compo = [ vector.Vector(v.compo[i:i+row]) for i in range(0, len(v.compo), row) ]
    return matrix.Matrix(row, col, compo)
    
def weaksim_worker((Ap, Bp)):
    try:
        A = unpack_matrix(Ap)
        if Bp is None: B = A
        else: B = unpack_matrix(Bp)
        lineq = MatrixLinearEquations(A.coeff_ring, A.row, A.column, parallel=False)
        lineq.weakSimilarity(A, B)
        return [ pack_matrix(M) for M in lineq.kernel() ]
    except KeyboardInterrupt:
        pass

def parallel_weaksim(As, Bs=None):
    if Bs is None: Bs = As
    return pool.pool().map_async(weaksim_worker, 
        [ (pack_matrix(A), pack_matrix(B)) for (A, B) in zip(As, Bs) ])

def intersect_worker(Lp):
    L = [ map(unpack_matrix, X) for X in Lp ]
    print "Intersecting %d solution spaces" % len(L)
    R = reduce(lambda x, y: intersect_solutions(x, y, False), L)
    return map(pack_matrix, R)

def intersect_solutions(S1, S2, parallel=True):
    if not S1 or not S2: return []
    field = S1[0].coeff_ring
    lineq = LinearEquations(field, len(S1) + len(S2), parallel=parallel)
    lineq.intersect(S1, S2)
    K = lineq.kernel()
    if K is None: return []
    Kp = [ vector.Vector(v.compo[0:len(S1)]) for v in K ]
    intersection = [ vector_to_matrix(k, S1) for k in Kp ]
    return [ X for X in intersection if X != matrix.zeroMatrix(X.row) ]

# Reduce the set of column vectors X to a basis
def basis_reduce(X):
    # We may select the first basis element arbitrarily
    B = matrix.Matrix(X.row, 1, [X[1]], X.coeff_ring)
    for col in range1(2, X.column):
        U = B.copy()
        U.extendColumn(X[col])
        if U.rank() > B.rank():
            B = U
    return B

# Convert the vector representation of a matrix over a basis B into
# its standard matrix representation M = sum_i (b_i * B_i)
def vector_to_matrix(v, B):
    M = matrix.Matrix(B[0].row, B[0].column, B[0].coeff_ring)
    for i in range1(len(v)):
        M += v[i] * B[i-1]
    return M

# Forward/back-substitution
def substitute(M, v, backward=False):
    x = vector.Vector([ M.coeff_ring.zero for i in range1(M.column) ])
    rows = list(range1(M.column))
    if backward: rows.reverse()
    for k in range(len(rows)):
        i = rows[k]
        y = v[i] - sum(x[j] * M[i,j] for j in rows[:k])
        x[i] = y * M[i,i].inverse()
    return x
    
# Decompose a matrix M into a linear combination of basis elements
# For basis {A_1, ..., A_n}, where M = sum_{k=1..n} a_k * A_k, return (a_k)
def decompose(M, (L, U, P)):
    v = flatten_matrix(M)
    return LUPsubstitute(v, (L, U, P))

def LUPsubstitute(v, (L, U, P)):
    # We want to find x s.t. Ax = v. We have that P * A = L * U.
    # Therefore we know L * U * x = P * v.
    # First compute w = P * v
    w = P * v
    # Next solve Ly = w by forward substitution
    y = substitute(L, w)
    #assert L * y == w
    # Then Ux = y by backward substitution
    x = substitute(U, y, backward=True)
    #assert U * x == y
    return x

def LUPsolve(A, b):
    LUP = A.LUPDecomposition()
    x = LUPsubstitute(b, LUP)
    assert A * x == b
    return x

# Given a matrix basis A_1 ... A_n, construct a matrix A = [A*1|...|A*n]
# where A*i is A_i flattened to a vector, and return (L, U, P, At) such that
# P * At * A = L * U, where P is a permutation matrix, and L and U are resp.
# upper and lower triangular.
def LUbasis(basis, field):
    dim = basis[0].row * basis[0].column
    A = matrix.Matrix(dim, len(basis),
            [ flatten_matrix(B) for B in basis ], field)
    (L, U, P) = A.LUPDecomposition()
    assert P * A == L * U
    return (L, U, P)

# Given a matrix X and list of indices l, return a matrix Y containing
# only the rows and columns of X whose indices are in l.    
def rc_eliminate(X, l):
    return X.subMatrix(l, l)

if __name__ == "__main__":
    A = matrix.FieldMatrix(2, 3, [5, 2, 1, 2, 3, 1], GF(7))
    b = ffvector([4, 1], GF(7))
    lineq = LinearEquations(A.coeff_ring, 3)
    lineq.addEquations(A, b)
    lineq.addEquation(ffvector([3, 1, 1], GF(7)), 2)
    assert (lineq.A*lineq.solution())[1] == ffvector(lineq.b.compo, A.coeff_ring)
    lineq.kernel()
    print A
    v = flatten_matrix(A)
    B = unflatten_matrix(v, A.row, A.column)
    assert A == B
    print
    print "Weak similarity"
    C = matrix.FieldMatrix(2, 2, [1, 0, 0, 1], GF2)
    D = matrix.FieldMatrix(3, 3, [1, 0, 1, 0, 1, 0, 1, 1, 0], GF2)
    lineq2 = MatrixLinearEquations(C.coeff_ring, 3, 2)
    lineq2.weakSimilarity(C, D)
    for X in lineq2.kernel():
        print X
        print
        assert X*C == D*X
    
