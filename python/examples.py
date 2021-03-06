from nzmath import matrix
import algebra
from util import GF2

def generate_algebras():
    algs = [(1, 1), (2, 1), (2, 2), (3, 1)]
    for x in algs:
        yield algebra_example(*x)

def generate_modules():
    mods = [(1, 1), (1, 2)]
    for x in mods:
        yield module_example(*x)

class AlgebraDescription:
    def __init__(self, field, stconsts, semisimple):
        self.stconsts = stconsts
        self.semisimple = semisimple
        self.field = field
        self.dim = stconsts[0].row
    def algebra(self):
        return algebra.Algebra(self.field, self.stconsts)
        
class ModuleDescription:
    def __init__(self, field, stconsts, algdesc):
        self.stconsts = stconsts
        self.algdesc = algdesc
        self.field = field
        self.dim = stconsts[0].row
    def module(self):
        A = self.algdesc.algebra()
        return algebra.Module(self.field, self.stconsts, A)

# Returns the nth associative algebra of dimension dim
def algebra_example(dim, n):
    # Algebras are defined over GF(2)
    field = GF2
    # All the algebras are unital
    stconsts = [matrix.unitMatrix(dim, field)]
    if dim == 1:
        if n == 1:
            # 1D, semisimple, 1 * 1 = 1
            return AlgebraDescription(field, stconsts, True)
    if dim == 2:
        if n == 1:
            # 2D, not semisimple, identity 1, aa = 1
            stconsts.append(matrix.Matrix(2, 2, [(0,1),(1,0)], field))
            return AlgebraDescription(field, stconsts, False)
        if n == 2:
            # 2D, not semisimple, identity 1, aa = 0
            stconsts.append(matrix.Matrix(2, 2, [(0,1),(0,0)], field))
            return AlgebraDescription(field, stconsts, False)
    if dim == 3:
        if n == 1:
            # 3D, semisimple, identity 1, aa = a, bb = b, ab = ba = 0
            stconsts.append(matrix.Matrix(3, 3, [(0,1,0),(0,1,0),(0,0,0)], field))
            stconsts.append(matrix.Matrix(3, 3, [(0,0,1),(0,0,0),(0,0,1)], field))
            return AlgebraDescription(field, stconsts, True)
    raise ValueError("No such example")

def module_example(dim, n):
    field = GF2
    stconsts = []
    if dim == 1:
        if n == 1:
            alg = algebra_example(1, 1)
            stconsts.append(matrix.Matrix(1, 1, [1], field))
            return ModuleDescription(field, stconsts, alg)
        if n == 2:
            alg = algebra_example(1, 1)
            stconsts.append(matrix.Matrix(1, 1, [0], field))
            return ModuleDescription(field, stconsts, alg)

