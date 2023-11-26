/*****************************************************************************************
                            qryomollapack.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QRYOMOLLAPACK_H
#define QRYOMOLLAPACK_H

#define FORTRAN_UNDERSCORE

#if defined(FORTRAN_UNDERSCORE)
 #define dsyev dsyev_
 #define dgelss dgelss_
 #define dgetrf dgetrf_
 #define dgetri dgetri_
 #define dgemv dgemv_
 #define dgemm dgemm_
#endif

extern "C" {
 
 void dgemv(char* trans,int* m, int* n,double* alpha,double* a,int* lda,double* x,int* incx,double* beta,double* y,int* incy);
 void dgetrf(int* m, int* n, double* a, int* lda,int* ipiv,int* info);
 void dgetri(int* n, double* a,int* lda,int* ipiv,double* work, int* lwork,int* info);
 void dsyev(char* jobz, char* uplo,int* n, double* a, int* lda, double* w, double* work, int* lwork, int* output);

 void dgelss(int* n, int* m, int* nhrs, double* a, int* lda, double* b, int* ldb, double* s, double* rcond,int* rank,double* work, int * lwork, int* info);
 void dgemm(char* transa, char* transb,int* m, int* n,int* k,double* alpha,double* a,int* lda,double* b, int* ldb,double* beta,double* c,int* ldc);
}

#endif //QRYOMOLLAPACK_H
