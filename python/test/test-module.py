from algebra import Module
from util import *
import examples
import pool

def test_reflexive():
    for Vd in mods:
        V = Vd.module()
        v = random_vector(V.dim, V.field)
        pi = V.reflexiveEndomorphism(v)
        assert (pi * v == v), "pi(v) != v"
        w = random_vector(V.algebra.dim, V.field)
        W = vector_to_matrix(w, V.stconsts)
        assert (W * pi == pi * W), "Not an endomorphism"

if __name__ == "__main__":
    mods = examples.generate_modules()
    pool.start_pool()
    test_reflexive()
    pool.stop_pool()
