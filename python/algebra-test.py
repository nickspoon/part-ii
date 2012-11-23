import algebra
from nzmath import finitefield, matrix

GF2 = finitefield.FinitePrimeField(2)
stconsts = []
stconsts.append(matrix.Matrix(2, 2, [1, 0, 0, 1], GF2))
stconsts.append(matrix.Matrix(2, 2, [0, 1, 1, 0], GF2))
L = algebra.Algebra(GF2, stconsts)
V = algebra.Module(GF2, L, stconsts)
print V.computeAnnihilators([1, 1])
print V.rankMax([0, 0])
