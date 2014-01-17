#ifndef _fixmatrixadditions_h_
#define _fixmatrixadditions_h_

#include "fixmatrix.h"

/*!
* \brief Calculates A*B*A'
*/
void mf16_mul_abat(mf16 *dest, const mf16 *a, const mf16 *b);

#endif