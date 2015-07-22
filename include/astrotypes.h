
#ifndef __ASTROTYPES_H__
#define __ASTROTYPES_H__

#define USE_LD

#ifdef USE_LD

typedef long double floatT;
#define COS_F cosl
#define SIN_F sinl
#define TAN_F tanl
#define SQRT_F sqrtl
#define POW_F powl
#define ATAN_F atanl
#define ASIN_F asinl
#define ACOS_F acosl
#define COSH_F coshl
#define SINH_F sinhl
#define TANH_F tanhl
#define STRTOD_F strtold


#else

typedef double floatT;
#define COS_F cos
#define SIN_F sin
#define TAN_F tan
#define SQRT_F sqrt
#define POW_F pow
#define ATAN_F atan
#define ASIN_F asin
#define ACOS_F acos
#define SINH_F sinh
#define COSH_F cosh
#define TANH_F tanh
#define STRTOD_F strtod

#endif

typedef int intT;

#endif

