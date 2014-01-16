#ifndef _fixmatrixadditions_h_
#define _fixmatrixadditions_h_

#include "fixmatrix.h"

/*!
* \brief Inverts a lower triangular square matrix
* \param[out] dest The destination matrix
* \param[in] matrix The lower triangular square matrix to invert
*
* Kudos: https://code.google.com/p/efficient-java-matrix-library
*/
void mf16_invert_lt(mf16 *dest, const mf16 *matrix);

#endif