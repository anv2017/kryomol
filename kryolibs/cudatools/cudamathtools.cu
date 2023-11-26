/*****************************************************************************************
                            cudamathtools.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/
#include <cuda.h>
#include <cuda_runtime.h>

#include "cudamathtools.h"

#ifndef M_PI
#define float M_PI 3.14159265358979323846;
#endif

#ifndef M_N
#define int M_N 256;
#endif


//Kernel functions
__global__ void VecAdd(float* a, float* b, float* c, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        c[i] = a[i]+b[i];
}

__global__ void VecSubs(float *a, float *b, float *c, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        c[i] = a[i]-b[i];

}

__global__ void ScMult(float w, float *a, float *b, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        b[i] = w*a[i];

}

__global__ void VecMult(float *a, float *b, float *c, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        c[i] = a[i]*b[i];

}

__global__ void Hamard(float *a, float *b, float *c, float *d, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] = a[i]*b[i]+c[i];
}

__global__ void SquareDiff(float w, float  *a, float *b, float *c, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        c[i] += w*(a[i]*a[i]-b[i]*b[i]);
}

__global__ void OrbitalS(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)));
}

__global__ void OrbitalPx(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(a[i]-cx);
}

__global__ void OrbitalPy(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(b[i]-cy);
}

__global__ void OrbitalPz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(c[i]-cz);
}

__global__ void OrbitalDxx(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(a[i]-cx)*(a[i]-cx);
}

__global__ void OrbitalDxy(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(a[i]-cx)*(b[i]-cy);
}

__global__ void OrbitalDxz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(a[i]-cx)*(c[i]-cz);
}

__global__ void OrbitalDyy(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(b[i]-cy)*(b[i]-cy);
}

__global__ void OrbitalDyz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(b[i]-cy)*(c[i]-cz);
}

__global__ void OrbitalDzz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(c[i]-cz)*(c[i]-cz);
}

__global__ void OrbitalFxxx(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(a[i]-cx)*(a[i]-cx)*(a[i]-cx);
}

__global__ void OrbitalFxxy(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(a[i]-cx)*(a[i]-cx)*(b[i]-cy);
}

__global__ void OrbitalFxxz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(a[i]-cx)*(a[i]-cx)*(c[i]-cz);
}

__global__ void OrbitalFxyy(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(a[i]-cx)*(b[i]-cy)*(b[i]-cy);
}

__global__ void OrbitalFxyz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(a[i]-cx)*(b[i]-cy)*(c[i]-cz);
}

__global__ void OrbitalFxzz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(a[i]-cx)*(c[i]-cz)*(c[i]-cz);
}

__global__ void OrbitalFyyy(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(b[i]-cy)*(b[i]-cy)*(b[i]-cy);
}

__global__ void OrbitalFyyz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(b[i]-cy)*(b[i]-cy)*(c[i]-cz);
}

__global__ void OrbitalFyzz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(b[i]-cy)*(c[i]-cz)*(c[i]-cz);
}

__global__ void OrbitalFzzz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(c[i]-cz)*(c[i]-cz)*(c[i]-cz);
}

__global__ void OrbitalDY0(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(0.25)*(sqrt(15/(M_PI)))*(2*(c[i]-cz)*(c[i]-cz)-(a[i]-cx)*(a[i]-cx)-(b[i]-cy)*(b[i]-cy));
}

__global__ void OrbitalDY1(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(-0.5)*(sqrt(15/(2*M_PI)))*(c[i]-cz)*(a[i]-cx);
}

__global__ void OrbitalDY2(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(0.5)*(sqrt(15/(2*M_PI)))*(c[i]-cz)*(a[i]-cx);
}

__global__ void OrbitalDY3(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(0.25)*(sqrt(15/(2*M_PI)))*((a[i]-cx)*(a[i]-cx)-(b[i]-cy)*(b[i]-cy));
}

__global__ void OrbitalDY4(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(0.25)*(sqrt(15/(2*M_PI)))*((a[i]-cx)*(a[i]-cx)-(b[i]-cy)*(b[i]-cy));
}

__global__ void OrbitalFY0(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(0.25)*(sqrt(7/(M_PI)))*(c[i]-cz)*(2*(c[i]-cz)*(c[i]-cz)-3*(a[i]-cx)*(a[i]-cx)-3*(b[i]-cy)*(b[i]-cy));
}

__global__ void OrbitalFY1(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(-0.125)*(sqrt(21/(M_PI)))*(a[i]-cx)*(4*(c[i]-cz)*(c[i]-cz)-(a[i]-cx)*(a[i]-cx)-(b[i]-cy)*(b[i]-cy));
}

__global__ void OrbitalFY2(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(0.125)*(sqrt(21/(M_PI)))*(a[i]-cx)*(4*(c[i]-cz)*(c[i]-cz)-(a[i]-cx)*(a[i]-cx)-(b[i]-cy)*(b[i]-cy));
}

__global__ void OrbitalFY3(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(0.25)*(sqrt(105/(2*M_PI)))*(c[i]-cz)*((a[i]-cx)*(a[i]-cx)-(b[i]-cy)*(b[i]-cy));
}

__global__ void OrbitalFY4(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(0.25)*(sqrt(105/(2*M_PI)))*(c[i]-cz)*((a[i]-cx)*(a[i]-cx)-(b[i]-cy)*(b[i]-cy));
}

__global__ void OrbitalFY5(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(-0.125)*(sqrt(35/M_PI))*(a[i]-cx)*((a[i]-cx)*(a[i]-cx)-(b[i]-cy)*(b[i]-cy));
}

__global__ void OrbitalFY6(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N)
        d[i] += xs*exp(-((a[i]-cx)*(a[i]-cx)+(b[i]-cy)*(b[i]-cy)+(c[i]-cz)*(c[i]-cz)))*(0.125)*(sqrt(35/M_PI))*(a[i]-cx)*((a[i]-cx)*(a[i]-cx)-(b[i]-cy)*(b[i]-cy));
}


//Cuda functions
void scudaadd(float *a, float *b, float *c, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    VecAdd<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,N);

    cudaMemcpy(c,dc,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
}

void scudasubs(float *a, float *b, float *c, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    VecSubs<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,N);

    cudaMemcpy(c,dc,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);

}

void scudascalmult(float w, float *a, float *b, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    ScMult<<<blocksPerGrid,threadsPerBlock>>>(w,da,db,N);

    cudaMemcpy(b,db,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);

}

void scudavecmult(float *a, float *b, float *c, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    VecMult<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,N);

    cudaMemcpy(c,dc,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);

}

void scudahamard(float *a, float *b, float *c, float *d, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    Hamard<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);

}

void scudasquarediff(float wa, float *a, float *b, float *c, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    SquareDiff<<<blocksPerGrid,threadsPerBlock>>>(wa,da,db,dc,N);

    cudaMemcpy(c,dc,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);

}

void scudaorbitals(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalS<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitalpx(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalPx<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitalpy(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalPy<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitalpz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalPz<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitaldxx(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalDxx<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitaldxy(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalDxy<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitaldxz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalDxz<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitaldyy(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalDyy<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitaldyz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalDyz<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitaldzz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalDzz<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitalfxxx(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalFxxx<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitalfxxy(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalFxxy<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitalfxxz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalFxxz<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitalfxyy(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalFxyy<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitalfxyz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalFxyz<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitalfxzz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalFxzz<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitalfyyy(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalFyyy<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitalfyyz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalFyyz<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitalfyzz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalFyzz<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitalfzzz(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalFzzz<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}


void scudaorbitaldy0(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalDY0<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitaldy1(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalDY1<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitaldy2(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalDY2<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitaldy3(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalDY3<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitaldy4(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalDY4<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitalfy0(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalFY0<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitalfy1(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalFY1<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitalfy2(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalFY2<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitalfy3(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalFY3<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitalfy4(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalFY4<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitalfy5(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalFY5<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}

void scudaorbitalfy6(float *a, float *b, float *c, float *d, float xs, float cx, float cy,float cz, int N)
{
    size_t size = N*sizeof(float);

    float *da;
    float *db;
    float *dc;
    float *dd;

    cudaMalloc(&da,size);
    cudaMalloc(&db,size);
    cudaMalloc(&dc,size);
    cudaMalloc(&dd,size);

    cudaMemcpy(da,a,size,cudaMemcpyHostToDevice);
    cudaMemcpy(db,b,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dc,c,size,cudaMemcpyHostToDevice);

    int threadsPerBlock;
    int blocksPerGrid;

    if (N>M_N)
    {
        threadsPerBlock = M_N;
        blocksPerGrid = ceil( (N+threadsPerBlock-1)/threadsPerBlock);
    }

    OrbitalFY6<<<blocksPerGrid,threadsPerBlock>>>(da,db,dc,dd,xs,cx,cy,cz,N);

    cudaMemcpy(d,dd,size,cudaMemcpyDeviceToHost);

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);
    cudaFree(dd);
}
