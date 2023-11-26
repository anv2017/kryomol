/** Subroutines from the ccmath 2.2.1 library*/

#include "mathtoolsapi.h"

namespace ccmath
{
MATHTOOLS_API void trnm(double *a,int m);
MATHTOOLS_API void housev(double *a,double *d,double *dp,int n);
MATHTOOLS_API int qrevec(double *ev,double *v,double *d,int m);
MATHTOOLS_API void eigen(double *a,double *ev,int n);
MATHTOOLS_API int minv(double *a,int n);
MATHTOOLS_API void vmul(double *vp,double *mat,double *v,int n);
MATHTOOLS_API void mmul(double *c,double *a,double *b,int n);
MATHTOOLS_API int solvps(double *a,double *b,int n);
MATHTOOLS_API int solv(double* a, double* b, int n);
MATHTOOLS_API void mattr(double *a,double *b,int m,int n);

MATHTOOLS_API void trnm(float *a,int m);
MATHTOOLS_API void housev(float *a,float *d,float *dp,int n);
MATHTOOLS_API int qrevec(float *ev,float *v,float *d,int m);
MATHTOOLS_API void eigen(float *a,float *ev,int n);
MATHTOOLS_API int minv(float *a,int n);
MATHTOOLS_API void vmul(float *vp,float *mat,float *v,int n);
MATHTOOLS_API void mmul(float *c,float *a,float *b,int n);
MATHTOOLS_API int solvps(float *a,float *b,int n);
MATHTOOLS_API int solv(float* a, float* b, int n);
MATHTOOLS_API void mattr(float *a,float *b,int m,int n);
MATHTOOLS_API double lsqsv(double *x,int *pr,double *var,double *d,double *b,double *v,
		int m,int n,double th);
MATHTOOLS_API void smgen(double *a,double *eval,double *evec,int n);
MATHTOOLS_API int sv2lsq(double *d,double *a,double *b,int m,double *v,int n);
MATHTOOLS_API int svdlsq(double *d,double *a,double *b,int m,double *v,int n);
MATHTOOLS_API void ldvmat(double *a,double *v,int n);
MATHTOOLS_API int qrbdbv(double *d,double *e,double *b,double *v,int n);
}
