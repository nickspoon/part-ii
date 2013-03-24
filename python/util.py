from nzmath import finitefield, vector, matrix, prime
from nzmath.rational import theIntegerRing
import numpy as np
import random
from math import ceil
#import functools

GF = lambda n: finitefield.FinitePrimeField(n)
GF2 = finitefield.FinitePrimeField(2)
next_prime_field = lambda n: GF(prime.nextPrime(n))

class DimensionError(vector.VectorSizeError, matrix.MatrixSizeError):
    pass

def matrixToSet(X):
    return set( X[i] for i in range1(X.column) )

def ffvector(compo, field):
    ffcmp = map(field.createElement, compo)
    return vector.Vector(ffcmp)
    
def ffmatrix(M, field):
    return M.map(field.createElement)

def range1(start, stop=None, step=None):
    if stop is None:
        return range(1, start+1)
    elif step is None:
        return range(start, stop+1)
    else:
        return range(start, stop+1, step)
        
def print_matrices(l):
    strs = [ str(A).split('\n') for A in l ]
    maxlen = max(len(x) for x in strs)
    stlen = len(strs[0][0])
    pad = "-" * stlen
    for x in strs:
        x += [ pad ] * (maxlen - len(x))
    for i in range(l[0].row):
        print "\t".join(x[i] for x in strs)

def random_element(field):
    p = field.getCharacteristic()
    if p == 0:
        return field.createElement(random.randint(0, 3))
    return field.createElement(random.randint(0, p))

def random_matrix(row, col, field=None):
    if field is None:
        field = col
        col = row
    R = matrix.Matrix(row, col, field)
    return R.map(lambda x: random_element(field))

def random_vector(dim, field):
    return vector.Vector([ random_element(field) for x in range1(dim) ])

def numpy_to_nzmath(arr, field):
    return matrix.Matrix(arr.shape[0], arr.shape[1], arr.tolist(), field)

def nzmath_to_numpy(M):
    Z = M.map(lambda x: x.getResidue())
    return np.array([x for x in Z.compo], np.int32)

# Given a matrix M and integer k, return M**k over Z using numpy
def numpy_matrix_pow(M, k):
    X = nzmath_to_numpy(M)
    Y = np.linalg.matrix_power(X, k)
    return numpy_to_nzmath(Y, theIntegerRing)
    
def pack_matrix(M):
    return (nzmath_to_numpy(M), M.coeff_ring.getCharacteristic())

def unpack_matrix((X, p)):
    return numpy_to_nzmath(X, GF(p))

def pack_vector(v):
    return ([ x.getResidue() for x in v.compo ], v[1].getModulus())

def unpack_vector((l,m)):
    field = GF(m)
    return vector.Vector([ field.createElement(x) for x in l ])

def numpy_save_matrices(matrices, fn):
    arrs = map(nzmath_to_numpy, matrices)
    arr = np.array(arrs)
    np.savez(fn, arr)

def numpy_load_matrices(fn, field):
    npz = np.load(fn)
    arr = npz['arr_0']
    matrices = []
    if len(arr.shape) == 2:
        arr = np.reshape(arr, (1, arr.shape[0], arr.shape[1]))
    for i in range(arr.shape[0]):
        matrices.append(numpy_to_nzmath(arr[i], field))
    npz.close()
    return matrices

def default_file_name(desc, field):
    return desc + "_" + str(field.getCharacteristic()) + ".npz"

# Given a list l, return a list of n equally-sized sublists of l
def nchunks(l, n):
    size = (len(l) + n - 1) / n
    return chunks(l, size)

# Given a list l, return a list of sublists of l of length size
def chunks(l, size):
    return (l[i:i+size] for i in range(0, len(l), size))

#def ki_wrap(f):
#    @functools.wraps(f)
#    def wrapped(*a, **k):
#        try:
#            f(*a, **k)
#        except KeyboardInterrupt:
#            pass
#    return wrapped
    
if __name__ == "__main__":
    print range1(10)
    print range1(2,10)
    print range1(2,10,2)
