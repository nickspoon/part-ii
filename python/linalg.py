from nzmath import matrix, vector
from util import *

# Class representing a system of linear equations.
# The solution is x where self.A * x = self.b
class LinearEquations:
    
    # A: m*n matrix, b: vector of length m
    def __init__(self, A, b):
        if A.row != len(b):
            raise ValueError("Matrix size mismatch")
        self.A = A
        self.b = b
        
    # X: k*n matrix; y: vector of length k
    def addEquations(self, X, y):
        if X.column != self.A.column:
            raise ValueError("Supplied matrix has the wrong number of columns")
        if X.row != len(y):
            raise ValueError("Matrix size mismatch")
        self.A.extendRow(X)
        self.b.compo.extend(y)
    
    def solve(self):
        return A.solve(b)

# Flatten a matrix column-wise into a vector
# e.g. [ [ a b ], [ c d ] ] => [ a c b d ]
def flatten_matrix(A):
    v = vector.Vector([])
    for i in range1(A.column):
        v.compo.extend(A[i])
    return v

if __name__ == "__main__":
    A = matrix.FieldMatrix(2, 3, [5, 2, 1, 2, 3, 1])
    b = vector.Vector([4, 1])
    lineq = LinearEquations(A, b)
    X = matrix.FieldMatrix(1, 3, [3, 1, 1])
    y = vector.Vector([2])
    lineq.addEquations(X, y)
    print lineq.solve()
    print A
    print flatten_matrix(A)
