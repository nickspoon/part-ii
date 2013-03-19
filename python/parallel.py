from nzmath import ring
from util import *
import types
import pool
import sys

def inject_parallel(M):
    M._cohensSimplify = types.MethodType(parallel_cohens, M)
    M.solve = types.MethodType(solve, M)

def parallel_cohens(self):
    """
    NS: Add some parallel processing to the cohensSimplify procedure
    _cohensSimplify is a common process used in image() and kernel()

    Return a tuple of modified matrix M, image data c and kernel data d.
    """
    print >>sys.stderr, "Parallel Cohens simplify %dx%d" % (self.row, self.column)
    M = self.copy()
    c = [0] * (M.row + 1)
    d = [-1] * (M.column + 1)
    Mpr = [ pack_vector(M.getRow(i)) for i in range1(M.row) ]
    for k in range(1, M.column + 1):
        for j in range(1, M.row + 1):
            if not c[j] and Mpr[j-1][0][k-1]:
                break
        else:           # not found j such that m(j, k)!=0 and c[j]==0
            d[k] = 0
            continue
        Mj = unpack_vector(Mpr[j-1])
        top = -ring.inverse(Mj[k])
        Mj[k] = -self.coeff_ring.one
        for s in range(k + 1, M.column + 1):
            Mj[s] = top * Mj[s]
        Mjp = pack_vector(Mj)
        Mpr[j-1] = Mjp
        work = [ (Mpr[i-1], Mjp, k)
                    for i in range(1, M.row + 1)
                    if i != j ]
        result = pool.pool().map(cohens_worker, work)
        i = 1
        for v in result:
            if i != j:
                Mpr[i-1] = v
            else:
                i += 1
                Mpr[i-1] = v
            i += 1
        c[j] = k
        d[k] = j
    for i in range(1, M.row + 1):
        M.setRow(i, unpack_vector(Mpr[i-1]))
    return (M, c, d)

def cohens_worker((Mip, Mjp, k)):
    Mi = unpack_vector(Mip)
    Mj = unpack_vector(Mjp)
    top = Mi[k]
    Mi[k] = Mi[1].getRing().zero
    for s in range1(k + 1, len(Mi)):
        Mi[s] = Mi[s] + top * Mj[s]
    return pack_vector(Mi)

def solve(self, B):  # modified Algorithm 2.3.4 of Cohen's book
    """
    Return solution X for self * X = B (B is a vector).
    This function returns tuple (V, M) below.
      V: one solution as vector
      M: kernel of self as list of basis vectors.
    If you want only one solution, use 'inverseImage'.

    Warning: B should not be a matrix instead of a vector
    """
    M_1 = self.copy()
    M_1.insertColumn(self.column + 1, B.compo)
    M_1._cohensSimplify = types.MethodType(parallel_cohens, M_1)
    V = M_1.kernel()
    ker = []
    flag = False
    if not V:
        raise NoInverseImage("no solution")
    n = V.row
    for j in range(1, V.column + 1):
        if not bool(V[n, j]): # self's kernel
            ker.append(vector.Vector([V[i, j] for i in range(1, n)]))
        elif not(flag):
            d = -ring.inverse(V[n, j])
            sol = vector.Vector([V[i, j] * d for i in range(1, n)])
            flag = True
    if not(flag):
        raise NoInverseImage("no solution")
    return sol, ker
