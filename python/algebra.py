from nzmath import finitefield, matrix, vector
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
                    try:
                        a = m.solve(b)
                        if a is not None:
                            return ann[z]
                    except matrix.NoInverseImage:
                        pass
        return None
        
    def reflexiveEndomorphism(self, v):
        # Calculate a projection pi which is a lambda-endomorphism of V, such 
        # that im(pi) = lambda*v and pi(v) = v by solving a system of linear
        # equations.
        
        # First guarantee that pi is a lambda-endomorphism
        m = matrix.Matrix(self.algebra.dim*self.dim**2, self.dim**2, self.field)
        for i in range(self.algebra.dim):
            for j in range(self.dim):
                for l in range(self.dim):
                    # Expand the equation indices
                    row = self.algebra.dim*self.dim*i + self.dim*j + l
                    for k in range(self.dim):
                        m[row, self.dim*k + j] += self.stconsts[i][k,l]
                        m[row, self.dim*l + k] -= self.stconsts[i][j,k]
        pi_space = m.kernel()
        pi = []
        for i in range(pi_space.column):
            pimat = matrix.Matrix(self.dim, self.dim, pi_space.getColumn(i).compo, \
                                    self.field)
            pi.append(pimat)
        return pi
