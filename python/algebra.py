from nzmath import matrix, vector
from linalg import MatrixLinearEquations
from util import *

class StructureConstantObject(object):
    def __init__(self, field, stconsts=None):
        self.stconsts = stconsts
        self.field = field
        if self.stconsts:
            self.dim = stconsts[0].column
        
    def getDimension(self):
        return self.dim

    # Get structure constant eta_i,j^k
    # Structure constants are stored as a list of matrices H_i.
    # > eta_i,j^k = H_i[k,j]
    # This representation allows for matrix multiplication to act on 
    # module-vectors as expected, i.e. the result of multiplying a module
    # element x with an algebra basis vector y is represented by V = H_y * V_x
    # where V_x is the vector representation of x.
    def getStruct(self, i, j, k):
        return self.stconsts[i-1][k,j]
        
    def vectorMultiply(self, a, b):
        M = matrix.Matrix(self.dim, self.dim, self.field)
        for i in range1(len(a)):
            M += a[i] * self.stconsts[i-1]
        return M*b
    
    # Return a spanning set S of Lv given by
    # S = { a_i * v | a_i in basis(L) }
    # N.B. for convenience this is not a set but a matrix, so that S[i] = a_i * v
    def spanningSet(self, v):
        return matrix.Matrix(self.dim, self.algebra.dim, [ s * v for s in self.stconsts ])

class Algebra(StructureConstantObject):
    pass

class Module(StructureConstantObject):
    
    def __init__(self, field, alg, stconsts=None):
        self.algebra = alg
        super(Module, self).__init__(field, stconsts)
        
    # This matrix represents the coefficients of a system of linear equations
    # of the form av = X, where a is an element of the algebra, and v and X
    # are elements of this module.
    # M[k,i] = sum[j](v[j]*stconsts[i,j,k])
    def multiplicativeMatrix(self, v):
        m = matrix.Matrix(self.dim, self.algebra.dim, self.field)
        for row in range1(m.row):
            for col in range1(m.column):
                m[row,col] = sum(v[j] * self.getStruct(col, j, row) for j in range1(self.dim))
        return m
    
    def computeAnnihilators(self, v, m=None):
        if m is None:
            m = self.multiplicativeMatrix(v)
        return m.kernel()
        
    def rankMax(self, v):
        # Find the first w in basis(V) such that Ann(V)w is not a subset of lambda*v
        # If v is of maximal rank, this procedure returns None, otherwise it returns w.
        # Save us computing the multiplicative matrix twice
        m = self.multiplicativeMatrix(v)
        ann = self.computeAnnihilators(v, m)
        if ann is not None:
            # Iterate through the basis of V v1...vn
            for j in range1(self.dim):
                # Iterate through the annihilators of v
                for z in range1(ann.column):
                    b = vector.Vector([ sum(ann[i,z] * self.getStruct(i, j, k) \
                            for i in range1(self.algebra.dim)) \
                            for k in range1(self.dim) ])
                    # ma = b has a solution iff av = Ann[z]w
                    try:
                        a = m.solve(b)
                    except matrix.NoInverseImage:
                        a = None
                    if a is None:
                        return ffvector([ 1 if x == j else 0 for x in range1(self.dim) ],
                                             self.field)
        return None
        
    def reflexiveEndomorphism(self, v):
        # Calculate a projection pi which is a lambda-endomorphism of V, such 
        # that im(pi) = lambda*v and pi(v) = v by solving a system of linear
        # equations.
        
        # First guarantee that pi is a lambda-endomorphism
        lineq = MatrixLinearEquations(self.field, self.dim, self.dim)
        for M in self.stconsts:
            lineq.weakSimilarity(M)
        
        # Add the condition pi(v) = v
        lineq.matrixEquation(v, v)
        solspace = lineq.solve()
        
        # Select a solution s.t. im(pi) = lambda*v
        return self.imageOf(v, solspace)
        
    def imageOf(self, v, (X, kernel)):
        # Given a matrix X and basis X_1, ..., X_n defining a space S, find
        # a projection pi in S such that im(pi) = lambda * v.
        
        # We require that for all v_i in basis(V), there is some x_i in lambda
        # such that pi(v_i) = x_i * v = sum_j(c_ij * a_j * v) for a_j in
        # basis(lambda).
        
        # From the precondition, we know that pi = X + sum_i(d_i * X_i) for 
        # some d_i.
        S = self.spanningSet(v)
        lineq = MatrixLinearEquations(self.field, self.dim, self.algebra.dim + 1)
        for row in range1(self.dim):
            for col in range1(self.dim):
                # Matrix of coefficients of c_ij, plus an extra column for d_i
                C = matrix.Matrix(self.dim, self.algebra.dim + 1, self.field)
                # sum_basis{(A[basis] * v)[row] * c[col,basis]}
                for basis in range1(self.algebra.dim):
                    C[col,basis] = S[row,basis]
                # - sum_i{(X_i)[row,col] * d_i}
                for i in range1(len(kernel)):
                    C[i,C.column] = -kernel[i-1][row,col]
                lineq.addCoefficientMatrix(C, X[row,col])
        soln = lineq.solution()
        d = soln[soln.column]
        Z = X.copy()
        for i in range1(len(kernel)):
            Z += d[i] * kernel[i-1]
        return Z
        
    def findMaxRank(self, initial):
        v = initial.copy()
        w = self.rankMax(v)
        while w is not None:
            pi = self.reflexiveEndomorphism(v)
            v += w - pi*w
            w = self.rankMax(v)
        return v
