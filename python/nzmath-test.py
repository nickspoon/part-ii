from nzmath import finitefield, matrix, vector
from nzmath.rational import theRationalField as R
from nzmath.poly import ring as polyring
from nzmath.poly.uniutil import polynomial
from util import *

f2 = finitefield.FinitePrimeField(2)
x = f2.createElement(0)
y = f2.createElement(1)
print x - y
print y + 1
print y.inverse()
print

f3 = finitefield.FinitePrimeField(3)
z = f3.createElement(2)
try:
    print y + z
except ValueError, e:
    print e.message
print

m = matrix.Matrix(2, 2, [1, 0, 1, 0], f2)
print m
print
v = matrix.Matrix(2, 1, [1, 1], f2)
print v
print
print m.inverseImage(v)
print

f2 = finitefield.FinitePrimeField(2)
P = polyring.PolynomialRing(f2)
x = P.createElement(polynomial([(1, 1)], f2))

def T(A, j):
    p = A.coeff_ring.getCharacteristic()
    phi = lambda A: (x * A + matrix.unitMatrix(A.row, A.coeff_ring)).determinant()
    return phi(A)[p**j]
    
A = matrix.Matrix(2, 2, [1, 0, 1, 0], f2)
print T(A, 1)
print

#A = matrix.Matrix(4, 4, [2, 1, 1, 1, 2, 1, 2, 1, 3, 2, 2, 1, 0, 1, 2, 1], R)
A = random_matrix(8, 3, f2)
print "Rank", A.rank()
try:
    (L, U, P) = A.LUPDecomposition()
    print "Result"
    B = L * U
    C = P * A
    print_matrices([L, U, P])
    print_matrices([A, C, B])
    assert C == B
except ZeroDivisionError:
    print "Div/zero"
    print A
