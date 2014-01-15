/*!
* \brief Microsoft Visual C keyword fix
*
* compile with "/FI msvcfix.h"
*/

#ifndef _MSVCFIX_H
#define _MSVCFIX_H

#ifdef _MSC_VER

#define inline __inline

#endif

#endif