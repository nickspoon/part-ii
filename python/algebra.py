from nzmath import finitefield, matrix

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
    
    def computeAnnihilators(self, v):
        m = matrix.Matrix(self.dim, self.algebra.dim, self.field)
        for row in range(0, m.row):
            for col in range(0, m.column):
                m[row,col] = sum(v[j] * self.stconsts[col][j,row] for j in range(0, self.dim))
        return m.kernel()
