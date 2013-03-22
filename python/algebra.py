from nzmath import matrix, vector
from linalg import *
from util import *
from itertools import product
import pool
import math
import sys
from parallel import inject_parallel

class StructureConstantObject(object):
    def __init__(self, field, stconsts):
        self.stconsts = stconsts
        self.field = field
        if self.stconsts:
            self.dim = stconsts[0].row
            self.unit = matrix.unitMatrix(self.dim, self.field)
        
    def getDimension(self):
        return self.dim

    # Get structure constant eta_i,j^k
    # Structure constants are stored as a list of matrices H_i.
    # > eta_i,j^k = H_i[k,j]
    # This representation allows for matrix multiplication to act on 
    # module-vectors as expected, i.e. the result of multiplying a module
    # element x with an algebra basis vector y is represented by V = H_y * V_x
    # where V_x is the vector representation of x.
    def getStruct(self, i, j, k):
        return self.stconsts[i-1][k,j]
        
    def toMatrix(self, a, B=None):
        if B is None: B = self.stconsts
        return vector_to_matrix(a, B)
        
    def vectorMultiply(self, a, b):
        return self.toMatrix(a)*b
        
    def basisListMatrix(self, l):
        return matrix.Matrix(self.dim, len(l), [ self.unit[i] for i in l ])
    
    # Return a spanning set S of Lv given by
    # S = { a_i * v | a_i in basis(L) }
    # N.B. for convenience this is not a set but a matrix, so that S[i] = a_i * v
    def spanningSet(self, v):
        return matrix.Matrix(self.dim, self.algebra.dim, [ s * v for s in self.stconsts ])
        
    # rk(v) is defined as rank(phi_v) where phi_v(x) = xv.
    # phi_v is equivalent to the matrix whose columns are a_i * v for a_i in basis(L)
    def rank(self, v):
        return self.spanningSet(v).rank()
        
    # Compute a basis B for this structure such that W + Z = B for 
    # Z < basis(self). Return [i] s.t. Z = <v_i> (i.e. self/span(W)).
    def quotientBasis(self, W):
        U = W.copy()
        I = self.unit
        bv = []
        for i in range1(self.dim):
            X = U.copy()
            X.extendColumn(I[i])
            inject_parallel(X)
            if X.rank() > U.rank():
                U = X
                bv.append(i)
        assert U.rank() == self.dim
        #Z = matrix.Matrix(self.dim, len(bv), [ I[i] for i in bv ], GF2)
        return bv
    
    # Computes structure constants for a new basis B of this structure
    # expressed as a matrix of vectors over the original basis
    def rebase(self, B1, B2=None):
        if B2 is None:
            B2 = B1
        # We express the multiplication of vectors in B1 by vectors in B2
        # in terms of a sum of the vectors in B2
        stconsts = []
        for i in range1(B1.column):
            M = matrix.Matrix(B2.column, B2.column, self.field)
            for j in range1(B2.column):
                v = self.vectorMultiply(B1[i], B2[j])
                M[j] = B2.inverseImage(v)[1]
            stconsts.append(M)
        assert all(B2 * stconsts[i-1][j] == self.vectorMultiply(B1[i], B2[j])
                    for i in range1(B1.column) for j in range1(B2.column))
        return stconsts

class Algebra(StructureConstantObject):

    @classmethod
    def fromBasis(cls, field, Lbasis):
        stconsts = struct_from_basis(field, Lbasis)
        return cls(field, stconsts)

    def __init__(self, field, stconsts):
        self.semisimple = False
        super(Algebra, self).__init__(field, stconsts)
    
    def radical(self):
        p = self.field.getCharacteristic()
        l = int(math.log(self.dim, p))
        # B = A union {1_n}
        # If we assume that the algebra is unital we don't need this
        B = self.stconsts + [matrix.unitMatrix(self.dim, self.field)]
        # I_-1 = A, represented by a matrix of basis vectors over A
        I = matrix.unitMatrix(len(B), self.field)
        I.deleteColumn(len(B))
        for i in range1(0, l):
            print >>sys.stderr, "Radical computation #%d" % i
            M = matrix.Matrix(len(B), I.column, self.field)
            inject_parallel(M)
            for k in range1(I.column):
                X = self.toMatrix(I[k], B)
                Xp = pack_matrix(X)
                work = [(p, i, Xp, pack_matrix(B[j-1])) for j in range1(len(B))]
                vec = pool.pool().map(compute_g_worker, work)
                M.setColumn(k, vector.Vector(vec))
                for j in range1(len(B)):
                    assert M[j,k] == compute_g(p, i, X*B[j-1])
            basis = M.kernel()
            if basis is None:
                self.semisimple = True
                return None
            # The solution is a set of column vectors sum_j(x_j * i_(i-1),j)
            # where i_k,1, ..., i_k,n = basis(I_k). Here we transform these 
            # back to vectors in the basis of B by multiplying by basis(I_k).
            zero = vector.Vector([ self.field.zero for i in range(I.row) ])
            I = matrix.Matrix(len(B), basis.column,
                    [ sum((basis[i,j] * I[i] for i in range1(I.column)), zero)
                        for j in range1(basis.column)])
        assert(all(x == self.field.zero for x in I.getRow(I.row)))
        return I.getBlock(1, 1, self.dim, I.column)
    
    # Given a basis for the quotient I, return the subalgebra L/I.
    def quotient(self, I):
        bv = self.quotientBasis(I)
        # Construct a basis B for self as vectors over the current basis
        B = matrix.Matrix(self.dim, len(bv), [ self.unit[i] for i in bv ])
        B.extendColumn(I)
        subconsts = self.rebase(B)[:len(bv)]
        for i in range(len(subconsts)):
            subconsts[i] = subconsts[i].getBlock(1, 1, len(bv))
        return Algebra(self.field, subconsts), self.basisListMatrix(bv)

class Module(StructureConstantObject):
    
    @classmethod
    def fromBasis(cls, field, Lbasis, Vbasis, alg):
        stconsts = struct_from_basis(field, Lbasis, Vbasis)
        return cls(field, stconsts, alg)
    
    def __init__(self, field, stconsts, alg):
        self.algebra = alg
        self.endo = None
        super(Module, self).__init__(field, stconsts)
        
    # This matrix represents the coefficients of a system of linear equations
    # of the form av = X, where a is an element of the algebra, and v and X
    # are elements of this module.
    # M[k,i] = sum[j](v[j]*stconsts[i,j,k])
    def multiplicativeMatrix(self, v):
        m = matrix.Matrix(self.dim, self.algebra.dim, self.field)
        for row in range1(m.row):
            for col in range1(m.column):
                m[row,col] = sum(v[j] * self.getStruct(col, j, row) for j in range1(self.dim))
        return m
    
    def computeAnnihilators(self, v, m=None):
        if m is None:
            m = self.multiplicativeMatrix(v)
        return m.kernel()
        
    def rankMax(self, v):
        # Find the first w in basis(V) such that Ann(V)w is not a subset of lambda*v
        # If v is of maximal rank, this procedure returns None, otherwise it returns w.
        
        # Save us computing the multiplicative matrix twice
        m = self.multiplicativeMatrix(v)
        ann = self.computeAnnihilators(v, m)
        if ann is not None:
            # Iterate through the basis of V v1...vn
            for j in range1(self.dim):
                # Iterate through the annihilators of v
                for z in range1(ann.column):
                    b = vector.Vector([ sum(ann[i,z] * self.getStruct(i, j, k) \
                            for i in range1(self.algebra.dim)) \
                            for k in range1(self.dim) ])
                    # ma = b has a solution iff av = Ann[z]w
                    try:
                        a = m.solve(b)
                    except matrix.NoInverseImage:
                        a = None
                    if a is None:
                        return self.unit[j]
        return None
    
    def endomorphismSpace(self):
        # Return a basis of the space of all endomorphisms of this module
        lineq = MatrixLinearEquations(self.field, self.dim, self.dim)
        for M in self.stconsts:
            lineq.weakSimilarity(M)
        self.endo = lineq.kernel()
        return self.endo
    
    def endomorphismSpace2(self):
        endspace = parallel_weaksim(self.stconsts)
        self.endo = ParallelIntersection(endspace.get(), packed=True).get()
        return self.endo
        
    def reflexiveEndomorphism(self, v):
        # Calculate a projection pi which is a lambda-endomorphism of V, such 
        # that im(pi) = lambda*v and pi(v) = v by solving a system of linear
        # equations.
        
        # First guarantee that pi is a lambda-endomorphism
        if self.endo is None:
            self.endomorphismSpace2()
        
        # Add the condition pi(v) = v
        lineq = LinearEquations(self.field, len(self.endo))
        lineq.matEqSpace(v, self.endo)
        (w, kernel) = lineq.solve()
        X = vector_to_matrix(w, self.endo)
        solspace = [ vector_to_matrix(k, self.endo)
                        for k in kernel ]
        
        # Select a solution s.t. im(pi) = lambda*v
        return self.imageOf(v, (X, solspace))
        
    def imageOf(self, v, (X, kernel)):
        # Given a matrix X and basis X_1, ..., X_n defining a space S, find
        # a projection pi in S such that im(pi) = lambda * v.
        
        # We require that for all v_i in basis(V), there is some x_i in lambda
        # such that pi(v_i) = x_i * v = sum_j(c_ij * a_j * v) for a_j in
        # basis(lambda).
        
        # From the precondition, we know that pi = X + sum_i(d_i * X_i) for 
        # some d_i.
        S = self.spanningSet(v)
        lineq = LinearEquations(self.field, self.dim * self.algebra.dim + len(kernel))
        for row in range1(self.dim):
            for col in range1(self.dim):
                # Matrix of coefficients of c_ij
                C = matrix.FieldMatrix(self.dim, self.algebra.dim, self.field)
                # sum_basis{(A[basis] * v)[row] * c[col,basis]}
                for basis in range1(self.algebra.dim):
                    C[col,basis] = S[row,basis]
                v = flatten_matrix(C)
                # - sum_i{(X_i)[row,col] * d_i}
                for i in range1(len(kernel)):
                    v.compo.append(-kernel[i-1][row,col])
                lineq.addEquation(v, X[row,col])
        #(soln, ker) = lineq.solve()
        soln = lineq.solution()
        d = vector.Vector(soln.compo[-len(kernel):])
        assert len(d) == len(kernel)
        Z = X.copy()
        for i in range1(len(kernel)):
            Z += d[i] * kernel[i-1]
        return Z
        
    def findGenerator(self):
        # If the algebra is not known to be semisimple, compute the radical
        if not self.algebra.semisimple:
            print >>sys.stderr, "Computing Rad(L)..."
            RadL = self.algebra.radical()
            if RadL is not None:
                print >>sys.stderr, RadL.column, "dimension(s)"
            else:
                print >>sys.stderr, "None"
            if RadL is not None:
                print >>sys.stderr, "Computing Rad(V) = Rad(L)V"
                RadV = self.radical(RadL)
                print >>sys.stderr, "Computing quotient algebra L/RadL...",
                Lbar, LB = self.algebra.quotient(RadL)
                print >>sys.stderr, Lbar.dim, "dimension(s)"
                # Lbar is semisimple since we have modded out its radical
                Lbar.semisimple = True
                print >>sys.stderr, "Computing quotient module V/RadV...",
                Vbar, VB = self.quotient(RadV, LB, Lbar)
                print Vbar.dim, "dimension(s)"
                print >>sys.stderr, "Finding a generator in V/RadV"
                x = Vbar.findGenerator()
                if x is not None:
                    return VB * x
                else:
                    return None
        # Semisimple case
        v = matrix.unitMatrix(self.dim, self.field)[1]
        print >>sys.stderr, "Have element v with rank", self.rank(v)
        w = self.rankMax(v)
        while w is not None:
            print >>sys.stderr, "Finding reflexive endomorphism pi"
            pi = self.reflexiveEndomorphism(v)
            assert pi * v == v
            v += w - pi*w
            print >>sys.stderr, "Have element v with rank", self.rank(v)
            w = self.rankMax(v)
        if self.rank(v) < self.dim: return None
        return v
        
    def radical(self, R=None):
        # RadV = (RadL)V = span({a_i*x_j | a_i in basis(RadL), x_j in basis(V)})
        if R is None: R = self.algebra.radical()
        if R is None: return None
        RV = None
        # a_i*x_j => matrix(a_i)*vector(x_j) = matrix(a_i)_j
        for col in range1(R.column):
            if RV is None:
                RV = self.toMatrix(R[col])
            else:
                RV.extendColumn(self.toMatrix(R[col]))
        return basis_reduce(RV)
    
    # Given the basis of an ideal I of V and the basis of the subalgebra over
    # which V/I is defined, return the submodule V/I.
    def quotient(self, I, Ba, alg):
        bv = self.quotientBasis(I)
        # Construct a basis B for self as vectors over the current basis
        B = matrix.Matrix(self.dim, len(bv),
            [ matrix.unitMatrix(self.dim, self.field)[i] for i in bv ])
        B.extendColumn(I)
        subconsts = self.rebase(Ba, B)
        for i in range(len(subconsts)):
            subconsts[i] = subconsts[i].getBlock(1, 1, len(bv))
        return Module(self.field, subconsts, alg), self.basisListMatrix(bv)

def compute_g(p, i, m):
    k = p**i
    # Compute the kth power of m as an integral matrix
    M = numpy_matrix_pow(m, k)
    # Compute f_i(M)
    f = M.trace()/k
    assert(f == int(f))
    # Obtain a result in F by taking modulo p
    return (int(f) % p)

def compute_g_worker((p, i, Xp, Bp)):
    X = unpack_matrix(Xp)
    B = unpack_matrix(Bp)
    return compute_g(p, i, X*B)
    
def decompose_worker((m1p, m2p, (Lp, Up, Pp))):
    (m1, m2, L, U, P) = map(unpack_matrix, (m1p, m2p, Lp, Up, Pp))
    result = m1 * m2
    v = decompose(result, (L, U, P))
    return pack_matrix(v.toMatrix(True))

# Computes the structure constant matrices for multiplication over bases
# basis1 = {A_1, ..., A_n}, basis2 = {V_1, ..., V_m} so that
# A_i * V_j = sum_k (eta_i,j^k * V_k)
def struct_from_basis(field, basis1, basis2=None):
    if basis2 is None: basis2 = basis1
    n = len(basis1)
    m = len(basis2)
    stconsts = [ matrix.Matrix(m, m, field) for i in range(n) ]
    LUb = LUbasis(basis2, field)
    packedLU = map(pack_matrix, LUb)
    packedb1 = map(pack_matrix, basis1)
    packedb2 = map(pack_matrix, basis2)
    vects = pool.pool().map(decompose_worker,
            [ (packedb1[i], packedb2[j], packedLU)
                for (i, j) in product(range(n), range(m)) ])
    unpacked = map(unpack_matrix, vects)
    for (i, j) in product(range(n), range(m)):
        stconsts[i].setColumn(j+1, unpacked[i * m + j][1])
        assert vector_to_matrix(stconsts[i][j+1], basis2) == basis1[i] * basis2[j]
    return stconsts
    
