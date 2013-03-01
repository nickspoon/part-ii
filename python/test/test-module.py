from algebra import Module
from util import *
import examples

def setUp():
    global mods
    mods = examples.generate_modules()

def test_reflexive():
    for Vd in mods:
        V = Vd.module()
        v = random_vector(V.dim, V.field)
        pi = reflexiveEndomorphism(v)
        assert (pi * v == v), "pi(v) != v"
        w = random_vector(V.algebra.dim, V.field)
