#include "typedefs.h"

class Algebra {
        Field K;
        std::vector<Matrix> basis;
        FieldArray3D structure;
        int dim;
        MDomain domain;
    public:
        Algebra(std::vector<Matrix> &basis, Field &K);
        Algebra(FieldArray3D &structure, int dim, Field &K);
        int dimension();
        FieldArray3D structureFromBasis();
};
