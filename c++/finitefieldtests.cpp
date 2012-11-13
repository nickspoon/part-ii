#include <iostream>
#include <string>
#include <linbox/field/modular.h>
#include <linbox/util/matrix-stream.h>
#include <linbox/blackbox/sparse.h>
#include <linbox/solutions/det.h>
#include <linbox/solutions/methods.h>
#include "typedefs.h"

using namespace std;
using namespace LinBox;

string MAT_STR = 
"2 2\n"
"1 2\n"
"3 4";

Matrix makeMatrix(Field F) {
    stringstream matstream;
    matstream << MAT_STR;
    MatrixStream<Field> ms(F, matstream);
    Matrix A(ms);
    return A;
}

int main(int argc, char** argv) {
    Field F(6);
    Matrix m = makeMatrix(F);
    MDomain D = MDomain(F);
    Field::Element d;
    det(d, m);
    cout << d << endl;
    D.add(m, m, m);
    det(d, m);
    cout << d << endl;
    return 0;
}
