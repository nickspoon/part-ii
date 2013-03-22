import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer, as_array
from util import *

ffpack = ctypes.cdll.LoadLibrary("./ffpack-interface.so")

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
    vec = nzmath_to_numpy(b.toMatrix())
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

if __name__ == "__main__":
    A = random_matrix(3, GF(3))
    b = random_vector(3, GF(3))
    K = kernel(A)
    assert K is None or all(A * K[i] == matrix.zeroMatrix(A.row, A.coeff_ring)[1]
                                for i in range1(K.column))
    v = solution(A, b)
    try:
        w = A.inverseImage(b)
    except matrix.NoInverseImage:
        w = None
    assert (v is None and w is None) or A * v == b
