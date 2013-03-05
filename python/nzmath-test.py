from nzmath import finitefield, matrix, vector
from nzmath.poly import ring as polyring
from nzmath.poly.uniutil import polynomial
import numpy as np

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

M = np.array([ f2.createElement(x) for x in [ 1, 0, 1, 0 ] ], dtype=object)
M = M.reshape(2,2)
print
print M + f2.one
