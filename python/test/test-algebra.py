from algebra import *
import examples

def setUp():
    global algs
    algs = examples.generate()

def test_radical():
    for A in algs:
        assert (A.algebra().radical() is None) == A.semisimple

def test_quotient():
    pass
