from nzmath import finitefield, matrix, vector

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
