from nzmath import finitefield, vector, matrix, prime
import numpy as np

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
    for i in range(l[0].row):
        print "\t".join(x[i] for x in strs)

def random_element(field):
    p = field.getCharacteristic()
    return field.createElement(random.randint(0, p))

def random_matrix(dim, field):
    R = matrix.Matrix(dim, dim, field)
    return R.map(lambda x: random_element(field))

def random_vector(dim, field):
    return vector.Vector([ random_element(field) for x in range1(dim) ])

def numpy_to_nzmath(arr, field):
    return matrix.Matrix(arr.shape[0], arr.shape[1], arr.tolist(), field)
        
if __name__ == "__main__":
    print range1(10)
    print range1(2,10)
    print range1(2,10,2)
