#include "fixmatrix.h"
#include "fixarray.h"

#include "fixmatrixadditions.h"

void mf16_mul_abat(mf16 *dest, const mf16 *a, const mf16 *b)
{
    mf16_mul(dest, a, b);
    mf16_mul_bt(dest, dest, a);
}
