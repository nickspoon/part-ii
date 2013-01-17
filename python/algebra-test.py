import algebra
from nzmath import finitefield, matrix, vector
from util import *

print "Module multiplication"
stconsts = []
stconsts.append(matrix.Matrix(2, 2, [1, 0, 0, 1], GF2))
stconsts.append(matrix.Matrix(2, 2, [0, 0, 1, 1], GF2))
L = algebra.Algebra(GF2)
V = algebra.Module(GF2, L, stconsts)
print V.vectorMultiply(ffvector([1, 0], GF2), ffvector([1, 1], GF2))
print

print "Module operations"
stconsts = []
stconsts.append(matrix.Matrix(3, 3, [1, 0, 0, 0, 1, 0, 0, 0, 1], GF2))
stconsts.append(matrix.Matrix(3, 3, [0, 0, 0, 1, 1, 0, 0, 0, 0], GF2))
stconsts.append(matrix.Matrix(3, 3, [0, 0, 0, 0, 0, 0, 1, 0, 1], GF2))
L = algebra.Algebra(GF2, stconsts)
stconsts = []
stconsts.append(matrix.Matrix(2, 2, [1, 0, 0, 1], GF2))
stconsts.append(matrix.Matrix(2, 2, [1, 0, 0, 0], GF2))
stconsts.append(matrix.Matrix(2, 2, [0, 0, 0, 1], GF2))
V = algebra.Module(GF2, L, stconsts)
v = ffvector([1, 0], GF2)
Ann = V.computeAnnihilators(v)
print "Annihilators"
print Ann
if Ann:
    for i in range1(Ann.column):
        assert V.vectorMultiply(Ann[i], v) == ffvector([0, 0], GF2)
print "w s.t. Ann(v)w not in Lv"
w = V.rankMax(v)
print w
print "Reflexive endomorphism"
pi = V.reflexiveEndomorphism(v)
print pi
assert pi * v == v
for x in ([0, 0], [1, 0], [0, 1], [1, 1]):
    for y in ([0, 0], [1, 0], [0, 1], [1, 1]):
        assert pi * V.vectorMultiply(ffvector(x, GF2), ffvector(y, GF2)) == \
                    V.vectorMultiply(ffvector(x, GF2), pi * ffvector(y, GF2)), \
                "pi * x * y != x * pi * y for x = " + str(x) + ", y = " + str(y)
print "An element of V which is of maximal rank"
print V.findMaxRank(v)

