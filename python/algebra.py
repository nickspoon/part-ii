from nzmath import finitefield, matrix, vector

class Algebra:

    def __init__(self, field, stconsts):
        self.stconsts = stconsts
        self.field = field
        self.dim = stconsts[0].column

class Module:
    
    def __init__(self, field, alg, stconsts):
        self.stconsts = stconsts
        self.algebra = alg
        self.field = field
        self.dim = stconsts[0].column

    def getDimension(self):
        return self.dim
        
    # This matrix represents the coefficients of a system of linear equations
    # of the form av = X, where a is an element of the algebra, and v and X
    # are elements of this module.
    # M[k,i] = sum[j](v[j]*stconsts[i,j,k])
    def multiplicativeMatrix(self, v):
        m = matrix.Matrix(self.dim, self.algebra.dim, self.field)
        for row in range(0, m.row):
            for col in range(0, m.column):
                m[row,col] = sum(v[j] * self.stconsts[col][j,row] for j in range(0, self.dim))
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
            for j in range(self.dim):
                # Iterate through the annihilators of v
                for z in range(ann.column):
                    b = vector.Vector([ sum(ann[z,i] * self.stconsts[i][j,k] \
                            for i in range(self.algebra.dim)) \
                            for k in range(self.dim) ])
                    try:
                        a = m.solve(b)
                        if a is not None:
                            return ann[z]
                    except matrix.NoInverseImage:
                        pass
        return None
