#include "fixmatrix.h"
#include "fixarray.h"

#include "fixmatrixadditions.h"

/*!
* \brief Inverts a lower triangular square matrix
* \param[out] dest The destination matrix
* \param[in] matrix The lower triangular square matrix to invert
*
* Kudos: https://code.google.com/p/efficient-java-matrix-library
*/
void mf16_invert_lt(mf16 *dest, const mf16 *matrix)
{
    int_fast8_t i, j, k;
    const uint_fast8_t n = matrix->rows;

    // If dest and input matrices alias, we have to use a temp matrix.
    mf16 tmp;
    fa16_unalias(dest, (void**)&dest, (void**)&matrix, &tmp, sizeof(tmp));

    dest->errors = dest->errors | matrix->errors;

    // TODO reorder these operations to avoid cache misses

    // inverts the lower triangular system and saves the result
    // in the upper triangle to minimize cache misses
    for (i = 0; i < n; ++i)
    {
        const fix16_t el_ii = matrix->data[i][i];
        for (j = 0; j <= i; ++j)
        {
            fix16_t sum = (i == j) ? fix16_one : 0;
            for (k = i - 1; k >= j; --k)
            {
                sum = fix16_sub(sum, fix16_mul(matrix->data[i][k], dest->data[j][k]));
            }
            dest->data[j][i] = fix16_div(sum, el_ii);
        }
    }
    // solve the system and handle the previous solution being in the upper triangle
    // takes advantage of symmetry
    for (i = n - 1; i >= 0; --i)
    {
        const fix16_t el_ii = matrix->data[i][i];
        for (j = 0; j <= i; ++j)
        {
            fix16_t sum = (i < j) ? 0 : dest->data[j][i];
            for (k = i + 1; k < n; ++k)
            {
                sum = fix16_sub(sum, fix16_mul(matrix->data[k][i], dest->data[j][k]));
            }
            dest->data[i][j] = dest->data[j][i] = fix16_div(sum, el_ii);
        }
    }
}