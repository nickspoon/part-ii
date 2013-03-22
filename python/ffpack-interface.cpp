#include <iostream>
#include "fflas-ffpack/field/modular-int32.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/utils/Matio.h"

using namespace std;
using namespace FFPACK;

typedef Modular<int32_t> Field;
typedef Field::Element Element;

/*
// Conversion is not needed here because the internal format is already
// an array of int32_t
Element* convert(const Field &F, const int *mat_in, int row, int col) {
    Element* M = new Element[row * col];
    int i = 0, j = 0;
    for (i = 0; i < row; i++) {
        for (j = 0; j < row; j++) {
            F.assign(*(M + i*row + j), *(mat_in + i*row + j));
        }
    }
    return M;
}*/

extern "C" int32_t* kernel(const int32_t *mat_in, const size_t row,
                         const size_t col, const int32_t p, size_t &kern_dim) {
    Field F(p);
    
    Element* M = (Element*)mat_in;
    Element* K;
    //write_field(F, cerr<<"M = "<<endl, M, row, col, col);
    size_t ldn;
    NullSpaceBasis(F, FFLAS::FflasRight, row, col, M, col, K, ldn, kern_dim);
    return K;
}

extern "C" int solve(const int32_t *mat_in, const int32_t* vect_in,
                        const size_t row, const size_t col, const int32_t p,
                        int32_t *result) {
    Field F(p);
    
    Element* M = (Element*)mat_in;
    Element* v = (Element*)vect_in;
    //write_field(F, cerr<<"M = "<<endl, M, row, col, col);
    //write_field(F, cerr<<"v = "<<endl, v, row, 1, 1);
    int info;
    fgesv(F, FFLAS::FflasLeft, row, col, 1, M, col, result, 1, v, 1, &info);
    //write_field(F, cerr<<"result = "<<endl, result, row, 1, 1);
    return info;
}
