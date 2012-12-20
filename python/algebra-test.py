import algebra
from nzmath import finitefield, matrix, vector
from util import *

GF2 = finitefield.FinitePrimeField(2)

print "Module multiplication"
stconsts = []
stconsts.append(matrix.Matrix(2, 2, [1, 0, 0, 1], GF2))
stconsts.append(matrix.Matrix(2, 2, [0, 0, 1, 1], GF2))
L = algebra.Algebra(GF2)
V = algebra.Module(GF2, L, stconsts)
print V.vectorMultiply(ffvector([1, 2], GF2), ffvector([1, 1], GF2))
print

print "Module operations"
stconsts = []
stconsts.append(matrix.Matrix(2, 2, [1, 0, 0, 1], GF2))
stconsts.append(matrix.Matrix(2, 2, [0, 1, 1, 0], GF2))
L = algebra.Algebra(GF2, stconsts)
V = algebra.Module(GF2, L, stconsts)
print V.computeAnnihilators(vector.Vector([1, 1]))
print V.rankMax(vector.Vector([1, 1]))
