from nzmath import finitefield, vector, matrix

GF = lambda n: finitefield.FinitePrimeField(n)
GF2 = finitefield.FinitePrimeField(2)

class DimensionError(vector.VectorSizeError, matrix.MatrixSizeError):
    pass

def matrixToSet(X):
    return set( X[i] for i in range1(X.column) )

def ffvector(compo, field):
    ffcmp = map(field.createElement, compo)
    return vector.Vector(ffcmp)

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
        
if __name__ == "__main__":
    print range1(10)
    print range1(2,10)
    print range1(2,10,2)
