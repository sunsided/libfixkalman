#ifndef _COMPILER_H_
#define _COMPILER_H_

/**
* \def RESTRICT Marks a restricted pointer 
*/
#ifdef __GNUC__
#define RESTRICT __restrict__
#else
#define RESTRICT
#endif

/**
* \def PURE Marks a function as pure, i.e. without global state
*/
#ifdef __GNUC__
#define PURE __attribute__ ((pure))
#else
#define PURE
#endif

/**
* \def CONST Marks a function as const, i.e. without global state and wihout using external memory
*/
#ifdef __GNUC__
#define CONST __attribute__ ((const))
#else
#define CONST
#endif

/**
* \def HOT Marks a function as a hot spot
*/
#ifdef __GNUC__
#define HOT __attribute__ ((hot))
#else
#define HOT
#endif

/**
* \def COLD Marks a function as a cold spot
*/
#ifdef __GNUC__
#define COLD __attribute__ ((cold))
#else
#define COLD
#endif

/**
* \def LEAF Marks a function as leaf
*/
#ifdef __GNUC__
#define LEAF __attribute__ ((leaf))
#else
#define LEAF
#endif

/**
* \def NONNULL Marks all function parameters as non-null
*/
#ifdef __GNUC__
#define NONNULL __attribute__ ((nonnull))
#else
#define NONNULL
#endif

/**
* \def RETURNS_NONNULL Marks the function return value as being non-null
*/
#ifdef __GNUC__
#define RETURNS_NONNULL __attribute__ ((returns_nonnull))
#else
#define RETURNS_NONNULL
#endif

/**
* \def INLINE Marks a function as to be inlined
*/
#ifdef __GNUC__
#define INLINE inline
#else
#ifdef _MSC_VER
#define INLINE __inline
#else
#define INLINE inline
#endif
#endif

/**
* \def EXTERN_INLINE Marks a function as to be inlined, but also externally defined
*/
#define EXTERN_INLINE extern INLINE

/**
* \def STATIC_INLINE Marks a function as to be inlined, but also statically defined
*/
#define STATIC_INLINE static INLINE

#endif