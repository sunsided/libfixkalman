#include <assert.h>
#include "fixmatrix.h"

#include "fixmatrixadditions.h"

/**
* \brief Tests matrix inversion using Cholesky decomposition
*/
void test_matrix_inverse()
{
    int result;

    // data buffer for the original and decomposed matrix
    mf16 d;
    d.columns = d.rows = 3;
    mf16_fill_diagonal(&d, fix16_one);
    d.data[0][1] = F16(0.5);
    d.data[1][0] = F16(0.5);
    
    // data buffer for the inverted matrix
    mf16 di;
    di.columns = di.rows = 3;
    mf16_fill(&di, 0);
    
    mf16 m;
    m.columns = m.rows = 3;

    // decompose matrix to lower triangular
    mf16_cholesky(&m, &d);

    // invert matrix using lower triangular
    mf16_invert_lt(&m, &m);

    // test the result
    fix16_t test = m.data[1][1];
    float testf = fix16_to_float(test);
    assert(testf >= 1.3);
}

void main() 
{
    test_matrix_inverse();
}