import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer, as_array
import nzmath.vector as vector
from os import path


def numpy_to_nzmath(arr, field):
    import nzmath.matrix as matrix
    return matrix.Matrix(arr.shape[0], arr.shape[1], arr.tolist(), field)

def nzmath_to_numpy(M):
    try:
        Z = M.map(lambda x: x.getResidue())
    except AttributeError:
        Z = M
    return np.array([x for x in Z.compo], np.int32)

def pack_matrix(M):
    return (nzmath_to_numpy(M), M.coeff_ring.getCharacteristic())

def unpack_matrix((X, p)):
    import nzmath.finitefield as finitefield
    GF = lambda n: finitefield.FinitePrimeField(n)
    return numpy_to_nzmath(X, GF(p))

ffpack_path = path.join(path.dirname(path.realpath(__file__)), "ffpack-interface.so")
ffpack = ctypes.cdll.LoadLibrary(ffpack_path)

def kernel(A):
    (arr, p) = pack_matrix(A)
    kern = ffpack.kernel
    kern.restype = ctypes.POINTER(ctypes.c_int)
    kern.argtypes = [ndpointer(ctypes.c_int), ctypes.c_size_t,
                    ctypes.c_size_t, ctypes.c_int,
                    ctypes.POINTER(ctypes.c_size_t)]
    kernel_size = ctypes.c_size_t(0)
    K_p = kern(arr, A.row, A.column, p,
            ctypes.byref(kernel_size))
    dim = kernel_size.value
    if dim == 0: return None
    k_arr = as_array(K_p, (A.column, dim))
    K = unpack_matrix((k_arr, p))
    return K

def solution(A, b):
    (arr, p) = pack_matrix(A)
    if isinstance(b, vector.Vector):
        vec = nzmath_to_numpy(b.toMatrix())
    else:
        vec = nzmath_to_numpy(b)
    solve = ffpack.solve
    solve.restype = ctypes.c_int
    solve.argtypes = [ndpointer(ctypes.c_int), ndpointer(ctypes.c_int),
                        ctypes.c_size_t, ctypes.c_size_t, ctypes.c_int,
                        ndpointer(ctypes.c_int)]
    vres = np.zeros((A.column, 1), dtype=np.int32)
    result = solve(arr, vec, A.row, A.column, p, vres)
    if result > 0: return None
    v = unpack_matrix((vres, p))
    return v

def LQUPDecomposition(A):
    (arr, p) = pack_matrix(A)
    lqup = ffpack.LQUP
    lqup.restype = ctypes.c_int
    lqup.argtypes = [ndpointer(ctypes.c_int), ctypes.c_size_t, ctypes.c_size_t,
        ctypes.c_int, ndpointer(ctypes.c_size_t), ndpointer(ctypes.c_size_t)]
    if ctypes.sizeof(ctypes.c_size_t) == 8:
        P = np.zeros(A.column, dtype=np.uint64)
        Q = np.zeros(A.row, dtype=np.uint64)
    else:
        P = np.zeros(A.column, dtype=np.uint32)
        Q = np.zeros(A.row, dtype=np.uint32)
    rank = lqup(arr, A.row, A.column, p, P, Q)
    LU = unpack_matrix((arr, p))
    return (LU, (P, Q), rank)

def LQUPSolve(LU, (P, Q), rank, b):
    (arr, p) = pack_matrix(LU)
    if isinstance(b, vector.Vector):
        vec = nzmath_to_numpy(b.toMatrix())
    else:
        vec = nzmath_to_numpy(b)
    solve = ffpack.LQUPsolve
    solve.restype = ctypes.c_int
    solve.argtypes = [ndpointer(ctypes.c_int), ndpointer(ctypes.c_int),
                        ctypes.c_size_t, ctypes.c_size_t,
                        ndpointer(ctypes.c_size_t), ndpointer(ctypes.c_size_t),
                        ctypes.c_int, ctypes.c_size_t,
                        ndpointer(ctypes.c_int)]
    vres = np.zeros((LU.column, 1), dtype=np.int32)
    result = solve(arr, vec, LU.row, LU.column, P, Q, p, rank, vres)
    if result > 0: return None
    v = unpack_matrix((vres, p))[1]
    return v
