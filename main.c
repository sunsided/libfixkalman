#include <assert.h>
#include "fixmatrix.h"

#if 0

/**
* \brief Tests matrix inversion using Cholesky decomposition
*/
void test_matrix_inverse()
{
    // the original matrix
    mf16 d = { 3, 3, 0, 0 };
    mf16_fill_diagonal(&d, fix16_one);
    d.data[0][1] = F16(0.5);
    d.data[1][0] = F16(0.5);

    // the decomposed and inverted matrix 
    mf16 m = { 3, 3, 0, 0 };
    mf16_fill(&m, 0);

    // decompose matrix to lower triangular
    mf16_cholesky(&m, &d);

    // invert matrix using lower triangular
    mf16_invert_lt(&m, &m);

    // test the result
    fix16_t test = m.data[1][1];
    float testf = fix16_to_float(test);
    assert(testf >= 1.3 && testf <= 1.4);
}

void main() 
{
    test_matrix_inverse();
}

#endif