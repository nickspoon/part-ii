from nzmath import matrix, vector
from nzmath import ring
from util import *

# Class representing a system of linear equations.
# The solution is x where self.A * x = self.b
class LinearEquations(object):
    
    # field: Field over which to solve
    # dim: number of variables
    def __init__(self, field, dim):
        self.field = field
        self.dim = dim
        self.A = None
        self.b = None
        
    # X: k*n matrix; y: vector of length k
    def addEquations(self, X, y):
        if X.column != self.dim:
            raise DimensionError("Supplied matrix has the wrong number of columns")
        if X.row != len(y):
            raise DimensionError("Matrix/vector size mismatch")
        if self.A is None:
            self.A = X.copy()
            self.b = y.copy()
        else:
            self.A.extendRow(X)
            self.b.compo.extend(y)
        
    # x: vector of length n; y: scalar
    def addEquation(self, x, y):
        X = matrix.Matrix(1, len(x), x.compo, self.field)
        Y = vector.Vector([y])
        self.addEquations(X, Y)
    
    # Returns a tuple (solution, kernel)
    def solve(self):
        return self.A.solve(self.b)
        
    # Returns a single solution
    def solution(self):
        return self.A.inverseImage(self.b)
    
    def kernel(self):
        return self.A.kernel()

# As above, but the variables form a matrix.
class MatrixLinearEquations(LinearEquations):
    def __init__(self, field, row, col):
        self.row = row
        self.col = col
        super(MatrixLinearEquations, self).__init__(field, row*col)
        
    def addCoefficientMatrix(self, A, n):
        self.addEquation(flatten_matrix(A), n)

    def weakSimilarity(self, A):
        # Find some X s.t. XA = AX, i.e. XA - AX = 0 where X, A are
        # matrices over a field.
        if A.column != A.row:
            raise DimensionError("Square matrix required for similarity")
        for i in range1(A.row):
            for j in range1(A.column):
                # Coefficients of X_ij, initially zero
                coeffs = matrix.FieldSquareMatrix(A.row, A.column, self.field)
                # XA: Multiply X-row by A-column
                coeffs.setRow(i, A.getColumn(j))
                # -AX: Multiply X-column by A-row
                coeffs.setColumn(j, -A.getRow(i))
                # XA - AX: handle row/column intersection
                coeffs[i,j] = A[j,j] - A[i,i]
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
        return self.solve()[0]
    
    def kernel(self):
        return self.solve()[1]

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

if __name__ == "__main__":
    A = matrix.FieldMatrix(2, 3, [5, 2, 1, 2, 3, 1])
    b = vector.Vector([4, 1])
    lineq = LinearEquations(A.coeff_ring, 3)
    lineq.addEquations(A, b)
    lineq.addEquation(vector.Vector([3, 1, 1]), 2)
    assert (lineq.A*lineq.solution())[1] == ffvector(lineq.b.compo, A.coeff_ring)
    print A
    v = flatten_matrix(A)
    B = unflatten_matrix(v, A.row, A.column)
    assert A == B
    print
    print "Weak similarity"
    C = matrix.FieldMatrix(2, 2, [1, 0, 0, 1], GF2)
    D = matrix.FieldMatrix(2, 2, [0, 1, 1, 0], GF2)
    lineq2 = MatrixLinearEquations(C.coeff_ring, 2, 2)
    lineq2.weakSimilarity(C)
    lineq2.weakSimilarity(D)
    empty = matrix.Matrix(2, 2, C.coeff_ring)
    for X in lineq2.kernel():
        print X
        assert (X*C - C*X) == empty
        assert (X*D - D*X) == empty
    
    print
    print "Ideal functions g_i"
    print compute_g(2, 1, D)
    
