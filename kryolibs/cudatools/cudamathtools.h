/*****************************************************************************************
                            cudamathtools.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef CUDAMATHTOOLS_H
#define CUDAMATHTOOLS_H

#include <cuda.h>
#include <cuda_runtime.h>
#include <cmath>

#ifdef WITH_CUDA

//Cuda functions
void scudaadd(float *a, float *b, float *c, int N);

void scudasubs(float *a, float *b, float *c, int N);

void scudascalmult(float w, float *a, float *b, int N);

void scudavecmult(float *a, float *b, float *c, int N);

void scudahamard(float *a, float *b, float *c, float *d, int N);

void scudasquarediff(float wa, float *a, float *b, float *c, int N);

void scudaorbitals(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitalpx(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitalpy(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitalpz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitaldxx(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitaldxy(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitaldxz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitaldyy(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitaldyz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitaldzz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitalfxxx(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitalfxxy(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitalfxxz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitalfxyy(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitalfxyz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitalfxzz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitalfyyy(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitalfyyz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitalfyzz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitalfzzz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitaldy0(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitaldy1(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitaldy2(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitaldy3(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitaldy4(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitalfy0(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitalfy1(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitalfy2(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitalfy3(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitalfy4(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitalfy5(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);

void scudaorbitalfy6(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N);


#endif


#endif // CUDAMATHTOOLS_H
