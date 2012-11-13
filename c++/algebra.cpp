#include "typedefs.h"
#include "algebra.h"

using namespace LinBox;
using namespace std;

Algebra::Algebra(vector<Matrix> &basis, Field &K)
    : basis(basis), K(K), domain(MDomain(K)), dim(basis.size()) {
    this->structure = this->structureFromBasis();
    //this->domain = MDomain(K);
    return;
}

Algebra::Algebra(FieldArray3D &structure, int dim, Field &K) 
    : structure(structure), K(K), dim(dim), domain(MDomain(K)) {
    return;
}

FieldArray3D Algebra::structureFromBasis() {
    FieldArray3D structure;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            Matrix C(this->K, dim, dim);
            this->domain.mul(C, basis[i], basis[j]);
            cout << C << endl;
        }
    }
    return structure;
}

int Algebra::dimension() {
    return this->dim;
}

Matrix matrixFromString(Field F, string str) {
    stringstream matstream;
    matstream << str;
    MatrixStream<Field> ms(F, matstream);
    Matrix A(ms);
    return A;
}

int main(int argc, char** argv) {
    Field F(3);
    FieldArray3D a;
    Algebra L(a, 0, F);
    vector<Matrix> b;
    string mat1 = 
    "2 2\n"
    "2 0\n"
    "0 2";
    b.push_back(matrixFromString(F, mat1));
    string mat2 = 
    "2 2\n"
    "0 1\n"
    "1 0";
    b.push_back(matrixFromString(F, mat2));
    Algebra M(b, F);
    cout << M.dimension() << endl;
}
