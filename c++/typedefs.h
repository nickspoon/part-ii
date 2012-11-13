#ifndef __TYPEDEFS_H__
#define __TYPEDEFS_H__
#include <linbox/field/modular.h>
#include <linbox/blackbox/sparse.h>
#include <linbox/matrix/blas-matrix.h>
typedef LinBox::Modular<short> Field;
typedef LinBox::BlasMatrix<Field> Matrix;
typedef LinBox::MatrixDomain<Field> MDomain;
typedef std::vector< std::vector< std::vector<Field::Element> > > FieldArray3D;
#endif
