/*****************************************************************************************
                            orbitalarray.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "orbitalarray.h"

using namespace kryomol;

void OrbitalArray::CalculateExponential(const float Nx, const float Ny, const float Nz, const float step, const float cx, const float cy, const float cz, const float alpha, const float xs) //const float &cx, const float &cy, const float &cz)
{
    assert ( this->m_rows>0 && this->m_columns>0 && this->m_files>0);
    #ifdef WITH_SIMD

    if ( (m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) )
    {
    __m128* ptr = (__m128*) m_data;
    __m128* v4 = (__m128*) _mm_malloc(4*sizeof(float),16);
    __m128* xs4 = (__m128*) _mm_malloc(4*sizeof(float),16);

    float a = -alpha*PC::AngstromToBohr*PC::AngstromToBohr;
    *xs4 = _mm_set_ps1(xs);

    #ifdef ROW_MAJOR
    size_t p=0;
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files/4;++k)
            {
                *v4 = _mm_set_ps(a*((x)*(x) + (y)*(y) + (z+3*step)*(z+3*step)),a*((x)*(x) + (y)*(y) + (z+2*step)*(z+2*step)),a*((x)*(x) + (y)*(y) + (z+step)*(z+step)),a*((x)*(x) + (y)*(y) + (z)*(z)));
                *v4 = exp_ps(*v4);
                ptr[p] = _mm_mul_ps(*xs4,*v4); //i*m_columns*m_files+j*m_files+i
                ++p;
                z += 4*step;
            }
            y += step;
        }
        x += step;
    }
    #else
    size_t p=0;
    float z = -Nz/2-cz;
    for(size_t k=0;k<m_files;++k)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            float x = -Nx/2-cx;
            for(size_t i=0;i<m_rows/4;++i)
            {
                *v4 = _mm_set_ps(a*((x+3*step)*(x+3*step) + (y)*(y) + (z)*(z)),a*((x+2*step)*(x+2*step) + (y)*(y) + (z)*(z)),a*((x+step)*(x+step) + (y)*(y) + (z)*(z)),a*((x)*(x) + (y)*(y) + (z)*(z)));
                *v4 = exp_ps(*v4);
                ptr[p] = _mm_mul_ps(*xs4,*v4); //k*m_rows*m_columns+j*m_rows+i
                ++p;
                x += 4*step;
            }
            y += step;
        }
        z += step;
    }
    #endif
    m_data = (float*) ptr;
    }
    else
    {
        float x = -Nx/2-cx;
        for(size_t i=0;i<m_rows;++i)
        {
            float y = -Ny/2-cy;
            for(size_t j=0;j<m_columns;++j)
            {
                float z = -Nz/2-cz;
                for(size_t k=0;k<m_files;++k)
                {
                    (*this)(i,j,k) += xs*(exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x)*(x) + (y)*(y) + (z)*(z))));
                    z += step;
                }
                y += step;
            }
            x += step;
        }
    }
#else
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) += xs*(exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x)*(x) + (y)*(y) + (z)*(z))));
                z += step;
            }
            y += step;
        }
        x += step;
    }
#endif
}


void OrbitalArray::CalculateOrbitalS( const OrbitalArray &exponential)
{
    assert ( this->m_rows == exponential.m_rows && this->m_columns == exponential.m_columns && this->m_files == exponential.m_files);
    #ifdef WITH_SIMD
        if ( (m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) )
        {
            __m128* ptr = (__m128*) (m_data);

            size_t j=0;
            for (size_t i=0; i<(m_rows*m_columns*m_files)/4; ++i)
            {
                ptr[i] = _mm_set_ps(exponential[j+3],exponential[j+2],exponential[j+1],exponential[j]);
                j+=4;
            }
            m_data = (float*) (ptr);
        }
        else
        {
            for(size_t i=0;i<m_rows;++i)
                for(size_t j=0;j<m_columns;++j)
                    for(size_t k=0;k<m_files;++k)
                        (*this)(i,j,k) = exponential(i,j,k);

        }
    #else
        for(size_t i=0;i<m_rows;++i)
        {
            for(size_t j=0;j<m_columns;++j)
            {
                for(size_t k=0;k<m_files;++k)
                {
                    (*this)(i,j,k) = exponential(i,j,k);
                }
            }
        }
        std::cout << std::endl;
    #endif
}


void OrbitalArray::CalculateOrbitalPx( const OrbitalArray &exponential, const float cx, const float Nx, const float step)
{
    assert ( this->m_rows == exponential.m_rows && this->m_columns == exponential.m_columns && this->m_files == exponential.m_files);

    #ifdef WITH_SIMD
        if  ((m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) )
        {
            __m128* ptr = (__m128*) (m_data);
            __m128* exp4 = (__m128*) (exponential.m_data);
            __m128* x4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
            #ifdef ROW_MAJOR
                size_t p=0;
                float x = -Nx/2-cx;
                for(size_t i=0;i<m_rows;++i)
                {
                    *x4 = _mm_set_ps1((x)*PC::AngstromToBohr);
                    for(size_t j=0;j<m_columns;++j)
                    {
                        for(size_t k=0;k<m_files/4;++k)
                        {
                            ptr[p] = _mm_mul_ps(*x4,exp4[p]);
                            ++p;
                        }
                    }
                    x+=step;
                }
            #else
                size_t p=0;
                for(size_t k=0;k<m_files;++k)
                {
                    for(size_t j=0;j<m_columns;++j)
                    {
                        float x = -Nx/2-cx;
                        for(size_t i=0;i<m_rows/4;++i)
                        {
                            *x4 = _mm_set_ps((x+3*step)*PC::AngstromToBohr,(x+2*step)*PC::AngstromToBohr,(x+step)*PC::AngstromToBohr,(x)*PC::AngstromToBohr);
                            ptr[p] = _mm_mul_ps(*x4,exp4[p]);
                            ++p;
                            x+=4*step;
                        }
                    }
                }
            #endif
        }
        else
        {
            float x = -Nx/2-cx;
            for(size_t i=0;i<m_rows;++i)
            {
                for(size_t j=0;j<m_columns;++j)
                {
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(x)*PC::AngstromToBohr;
                    }
                }
                x+=step;
            }
        }
    #else
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        for(size_t j=0;j<m_columns;++j)
        {
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*(x)*PC::AngstromToBohr;
            }
        }
        x+=step;
    }
#endif
}


void OrbitalArray::CalculateOrbitalPy( const OrbitalArray &exponential, const float cy, const float Ny, const float step)
{
    assert ( this->m_rows == exponential.m_rows && this->m_columns == exponential.m_columns && this->m_files == exponential.m_files);

    #ifdef WITH_SIMD
        if  ((m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) )
        {
            __m128* ptr = (__m128*) (m_data);
            __m128* exp4 = (__m128*) (exponential.m_data);
            __m128* y4 = (__m128*) (_mm_malloc(4*sizeof(float),16));

            #ifdef ROW_MAJOR
                size_t p=0;
                for(size_t i=0;i<m_rows;++i)
                {
                    float y = -Ny/2-cy;
                    for(size_t j=0;j<m_columns;++j)
                    {
                        *y4 = _mm_set_ps1((y)*PC::AngstromToBohr);
                        for(size_t k=0;k<m_files/4;++k)
                        {
                            ptr[p] = _mm_mul_ps(*y4,exp4[p]);
                            ++p;
                        }
                        y+=step;
                    }
                }
            #else
                size_t p=0;
                for(size_t k=0;k<m_files;++k)
                {
                    float y = -Ny/2-cy;
                    for(size_t j=0;j<m_columns;++j)
                    {
                        *y4 = _mm_set_ps1((y)*PC::AngstromToBohr);
                        for(size_t i=0;i<m_rows/4;++i)
                        {
                            ptr[p] = _mm_mul_ps(*y4,exp4[p]);
                            ++p;
                        }
                        y+=step;
                    }
                }
            #endif
        }
        else
        {
            for(size_t i=0;i<m_rows;++i)
            {
                float y = -Ny/2-cy;
                for(size_t j=0;j<m_columns;++j)
                {
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(y)*PC::AngstromToBohr;
                    }
                    y+=step;
                }
            }
        }
    #else
        for(size_t i=0;i<m_rows;++i)
        {
            float y = -Ny/2-cy;
            for(size_t j=0;j<m_columns;++j)
            {
                for(size_t k=0;k<m_files;++k)
                {
                    (*this)(i,j,k) = exponential(i,j,k)*(y)*PC::AngstromToBohr;
                }
                y+=step;
            }
        }
    #endif

}


void OrbitalArray::CalculateOrbitalPz( const OrbitalArray &exponential, const float cz, const float Nz, const float step)
{
    assert ( this->m_rows == exponential.m_rows && this->m_columns == exponential.m_columns && this->m_files == exponential.m_files);

    #ifdef WITH_SIMD
        if  ((m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) )
        {
            __m128* ptr = (__m128*) (m_data);
            __m128* exp4 = (__m128*) (exponential.m_data);
            __m128* z4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
            #ifdef ROW_MAJOR
                size_t p=0;
                for(size_t i=0;i<m_rows;++i)
                {
                    for(size_t j=0;j<m_columns;++j)
                    {
                        float z = -Nz/2-cz;
                        for(size_t k=0;k<m_files/4;++k)
                        {
                            *z4 = _mm_set_ps((z+3*step)*PC::AngstromToBohr,(z+2*step)*PC::AngstromToBohr,(z+step)*PC::AngstromToBohr,(z)*PC::AngstromToBohr);
                            ptr[p] = _mm_mul_ps(*z4,exp4[p]);
                            ++p;
                            z+=4*step;
                        }
                    }
                }
            #else
                size_t p=0;
                float z = -Nz/2-cz;
                for(size_t k=0;k<m_files;++k)
                {
                    *z4 = _mm_set_ps1((z)*PC::AngstromToBohr);
                    for(size_t j=0;j<m_columns;++j)
                    {
                        for(size_t i=0;i<m_rows/4;++i)
                        {
                            ptr[p] = _mm_mul_ps(*z4,exp4[p]);
                            ++p;
                        }
                    }
                    z+=step;
                }
            #endif
        }
        else
        {
            for(size_t i=0;i<m_rows;++i)
            {
                for(size_t j=0;j<m_columns;++j)
                {
                    float z = -Nz/2-cz;
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(z)*PC::AngstromToBohr;
                        z+=step;
                    }
                }
            }
        }
    #else
        for(size_t i=0;i<m_rows;++i)
        {
            for(size_t j=0;j<m_columns;++j)
            {
                float z = -Nz/2-cz;
                for(size_t k=0;k<m_files;++k)
                {
                    (*this)(i,j,k) = exponential(i,j,k)*(z)*PC::AngstromToBohr;
                    z+=step;
                }
            }
        }
    #endif
}


void OrbitalArray::CalculateOrbitalDxx( const OrbitalArray &exponential, const float cx, const float Nx, const float step)
{
#ifdef WITH_SIMD
    if  ((m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) )
    {
        __m128* ptr = (__m128*) (m_data);
        __m128* exp4 = (__m128*) (exponential.m_data);
        __m128* x4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        #ifdef ROW_MAJOR
            size_t p=0;
            float x = -Nx/2-cx;
            for(size_t i=0;i<m_rows;++i)
            {
                *x4 = _mm_set_ps1((x)*PC::AngstromToBohr);
                *x4 = _mm_mul_ps(*x4,*x4);
                for(size_t j=0;j<m_columns;++j)
                {
                    for(size_t k=0;k<m_files/4;++k)
                    {
                        ptr[p] = _mm_mul_ps(*x4,exp4[p]);
                        ++p;
                    }
                }
                x+=step;
            }
        #else
            size_t p=0;
            for(size_t k=0;k<m_files;++k)
            {
                for(size_t j=0;j<m_columns;++j)
                {
                    float x = -Nx/2-cx;
                    for(size_t i=0;i<m_rows/4;++i)
                    {
                        *x4 = _mm_set_ps((x+3*step)*PC::AngstromToBohr,(x+2*step)*PC::AngstromToBohr,(x+step)*PC::AngstromToBohr,(x)*PC::AngstromToBohr);
                        *x4 = _mm_mul_ps(*x4,*x4);
                        ptr[p] = _mm_mul_ps(*x4,exp4[p]);
                        ++p;
                        x+=4*step;
                    }
                }
            }
        #endif
    }
    else
    {
        float x = -Nx/2-cx;
        for(size_t i=0;i<m_rows;++i)
        {
            for(size_t j=0;j<m_columns;++j)
            {
                for(size_t k=0;k<m_files;++k)
                {
                    (*this)(i,j,k) = exponential(i,j,k)*(x)*PC::AngstromToBohr*(x)*PC::AngstromToBohr;
                }
            }
            x+=step;
        }
    }
#else
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        for(size_t j=0;j<m_columns;++j)
        {
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*(x)*PC::AngstromToBohr*(x)*PC::AngstromToBohr;
            }
        }
        x+=step;
    }
#endif
}


void OrbitalArray::CalculateOrbitalDyy( const OrbitalArray &exponential, const float cy, const float Ny, const float step)
{
#ifdef WITH_SIMD
    if  ((m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) )
    {
        __m128* ptr = (__m128*) (m_data);
        __m128* exp4 = (__m128*) (exponential.m_data);
        __m128* y4 = (__m128*) (_mm_malloc(4*sizeof(float),16));

        #ifdef ROW_MAJOR
            size_t p=0;
            for(size_t i=0;i<m_rows;++i)
            {
                float y = -Ny/2-cy;
                for(size_t j=0;j<m_columns;++j)
                {
                    *y4 = _mm_set_ps1((y)*PC::AngstromToBohr);
                    *y4 = _mm_mul_ps(*y4,*y4);
                    for(size_t k=0;k<m_files/4;++k)
                    {
                        ptr[p] = _mm_mul_ps(*y4,exp4[p]);
                        ++p;
                    }
                    y+=step;
                }
            }
        #else
            size_t p=0;
            for(size_t k=0;k<m_files;++k)
            {
                float y = -Ny/2-cy;
                for(size_t j=0;j<m_columns;++j)
                {
                    *y4 = _mm_set_ps1((y)*PC::AngstromToBohr);
                    *y4 = _mm_mul_ps(*y4,*y4);
                    for(size_t i=0;i<m_rows/4;++i)
                    {
                        ptr[p] = _mm_mul_ps(*y4,exp4[p]);
                        ++p;
                    }
                    y+=step;
                }
            }
        #endif
    }
    else
    {
        for(size_t i=0;i<m_rows;++i)
        {
            float y = -Ny/2-cy;
            for(size_t j=0;j<m_columns;++j)
            {
                for(size_t k=0;k<m_files;++k)
                {
                    (*this)(i,j,k) = exponential(i,j,k)*(y)*PC::AngstromToBohr*(y)*PC::AngstromToBohr;
                }
                y+=step;
            }
        }
    }
#else
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*(y)*PC::AngstromToBohr*(y)*PC::AngstromToBohr;
            }
            y+=step;
        }
    }
#endif
}


void OrbitalArray::CalculateOrbitalDzz( const OrbitalArray &exponential, const float cz, const float Nz, const float step)
{
#ifdef WITH_SIMD
    if  ((m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) )
    {
        __m128* ptr = (__m128*) (m_data);
        __m128* exp4 = (__m128*) (exponential.m_data);
        __m128* z4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        #ifdef ROW_MAJOR
            size_t p=0;
            for(size_t i=0;i<m_rows;++i)
            {
                for(size_t j=0;j<m_columns;++j)
                {
                    float z = -Nz/2-cz;
                    for(size_t k=0;k<m_files/4;++k)
                    {
                        *z4 = _mm_set_ps((z+3*step)*PC::AngstromToBohr,(z+2*step)*PC::AngstromToBohr,(z+step)*PC::AngstromToBohr,(z)*PC::AngstromToBohr);
                        *z4 = _mm_mul_ps(*z4,*z4);
                        ptr[p] = _mm_mul_ps(*z4,exp4[p]);
                        ++p;
                        z+=4*step;
                    }
                }
            }
        #else
            size_t p=0;
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                *z4 = _mm_set_ps1((z)*PC::AngstromToBohr);
                *z4 = _mm_mul_ps(*z4,*z4);
                for(size_t j=0;j<m_columns;++j)
                {
                    for(size_t i=0;i<m_rows/4;++i)
                    {
                        ptr[p] = _mm_mul_ps(*z4,exp4[p]);
                        ++p;
                    }
                }
                z+=step;
            }
        #endif
    }
    else
    {
        for(size_t i=0;i<m_rows;++i)
        {
            for(size_t j=0;j<m_columns;++j)
            {
                float z = -Nz/2-cz;
                for(size_t k=0;k<m_files;++k)
                {
                    (*this)(i,j,k) = exponential(i,j,k)*(z)*PC::AngstromToBohr*(z)*PC::AngstromToBohr;
                    z+=step;
                }
            }
        }
    }
#else
    for(size_t i=0;i<m_rows;++i)
    {
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*(z)*PC::AngstromToBohr*(z)*PC::AngstromToBohr;
                z+=step;
            }
        }
    }
#endif
}


void OrbitalArray::CalculateOrbitalDxy( const OrbitalArray &exponential, const float cx, const float cy, const float Nx, const float Ny, const float step)
{
#ifdef WITH_SIMD
    if  ((m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) )
    {
        __m128* ptr = (__m128*) (m_data);
        __m128* exp4 = (__m128*) (exponential.m_data);
        __m128* x4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        __m128* y4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        #ifdef ROW_MAJOR
            size_t p=0;
            float x = -Nx/2-cx;
            for(size_t i=0;i<m_rows;++i)
            {
                *x4 = _mm_set_ps1((x)*PC::AngstromToBohr);
                float y = -Ny/2-cy;
                for(size_t j=0;j<m_columns;++j)
                {
                    *y4 = _mm_set_ps1((y)*PC::AngstromToBohr);
                    *y4 = _mm_mul_ps(*x4,*y4);
                    for(size_t k=0;k<m_files/4;++k)
                    {
                        ptr[p] = _mm_mul_ps(*y4,exp4[p]);
                        ++p;
                    }
                    y+=step;
                }
                x+=step;
            }
        #else
            size_t p=0;
            for(size_t k=0;k<m_files;++k)
            {
                float y = -Ny/2-cy;
                for(size_t j=0;j<m_columns;++j)
                {
                    *y4 = _mm_set_ps1((y)*PC::AngstromToBohr);
                    float x = -Nx/2-cx;
                    for(size_t i=0;i<m_rows/4;++i)
                    {
                        *x4 = _mm_set_ps((x+3*step)*PC::AngstromToBohr,(x+2*step)*PC::AngstromToBohr,(x+step)*PC::AngstromToBohr,(x)*PC::AngstromToBohr);
                        *x4 = _mm_mul_ps(*x4,*y4);
                        ptr[p] = _mm_mul_ps(*x4,exp4[p]);
                        ++p;
                        x+=4*step;
                    }
                    y+=step;
                }
            }
        #endif
    }
    else
    {
        float x = -Nx/2-cx;
        for(size_t i=0;i<m_rows;++i)
        {
            float y = -Ny/2-cy;
            for(size_t j=0;j<m_columns;++j)
            {
                for(size_t k=0;k<m_files;++k)
                {
                    (*this)(i,j,k) = exponential(i,j,k)*(x)*PC::AngstromToBohr*(y)*PC::AngstromToBohr;
                }
                y+=step;
            }
            x+=step;
        }
    }
#else
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*(x)*PC::AngstromToBohr*(y)*PC::AngstromToBohr;
            }
            y+=step;
        }
        x+=step;
    }
#endif
}


void OrbitalArray::CalculateOrbitalDxz( const OrbitalArray &exponential, const float cx, const float cz, const float Nx, const float Nz, const float step)
{
#ifdef WITH_SIMD
    if  ((m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) )
    {
        __m128* ptr = (__m128*) (m_data);
        __m128* exp4 = (__m128*) (exponential.m_data);
        __m128* x4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        __m128* z4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        #ifdef ROW_MAJOR
            size_t p=0;
            float x = -Nx/2-cx;
            for(size_t i=0;i<m_rows;++i)
            {
                *x4 = _mm_set_ps1((x)*PC::AngstromToBohr);
                for(size_t j=0;j<m_columns;++j)
                {
                    float z = -Nz/2-cz;
                    for(size_t k=0;k<m_files/4;++k)
                    {
                        *z4 = _mm_set_ps((z+3*step)*PC::AngstromToBohr,(z+2*step)*PC::AngstromToBohr,(z+step)*PC::AngstromToBohr,(z)*PC::AngstromToBohr);
                        *z4 = _mm_mul_ps(*x4,*z4);
                        ptr[p] = _mm_mul_ps(*z4,exp4[p]);
                        ++p;
                        z+=4*step;
                    }
                }
                x+=step;
            }
        #else
            size_t p=0;
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                *z4 = _mm_set_ps1((z)*PC::AngstromToBohr);
                for(size_t j=0;j<m_columns;++j)
                {
                    float x = -Nx/2-cx;
                    for(size_t i=0;i<m_rows/4;++i)
                    {
                        *x4 = _mm_set_ps((x+3*step)*PC::AngstromToBohr,(x+2*step)*PC::AngstromToBohr,(x+step)*PC::AngstromToBohr,(x)*PC::AngstromToBohr);
                        *x4 = _mm_mul_ps(*x4,*z4);
                        ptr[p] = _mm_mul_ps(*x4,exp4[p]);
                        ++p;
                        x+=4*step;
                    }
                }
                z+=step;
            }
        #endif
    }
    else
    {
        float x = -Nx/2-cx;
        for(size_t i=0;i<m_rows;++i)
        {
            for(size_t j=0;j<m_columns;++j)
            {
                float z = -Nz/2-cz;
                for(size_t k=0;k<m_files;++k)
                {
                    (*this)(i,j,k) = exponential(i,j,k)*(x)*PC::AngstromToBohr*(z)*PC::AngstromToBohr;
                    z+=step;
                }
            }
            x+=step;
        }
    }
#else
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*(x)*PC::AngstromToBohr*(z)*PC::AngstromToBohr;
                z+=step;
            }
        }
        x+=step;
    }
#endif
}


void OrbitalArray::CalculateOrbitalDyz( const OrbitalArray &exponential, const float cy, const float cz, const float Ny, const float Nz, const float step)
{
#ifdef WITH_SIMD
    if  ((m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) )
    {
        __m128* ptr = (__m128*) (m_data);
        __m128* exp4 = (__m128*) (exponential.m_data);
        __m128* y4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        __m128* z4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        #ifdef ROW_MAJOR
            size_t p=0;
            for(size_t i=0;i<m_rows;++i)
            {
                float y = -Ny/2-cy;
                for(size_t j=0;j<m_columns;++j)
                {
                    *y4 = _mm_set_ps1((y)*PC::AngstromToBohr);
                    float z = -Nz/2-cz;
                    for(size_t k=0;k<m_files/4;++k)
                    {
                        *z4 = _mm_set_ps((z+3*step)*PC::AngstromToBohr,(z+2*step)*PC::AngstromToBohr,(z+step)*PC::AngstromToBohr,(z)*PC::AngstromToBohr);
                        *z4 = _mm_mul_ps(*z4,*y4);
                        ptr[p] = _mm_mul_ps(*z4,exp4[p]);
                        ++p;
                        z+=4*step;
                    }
                    y+=step;
                }
            }
        #else
            size_t p=0;
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                *z4 = _mm_set_ps1((z)*PC::AngstromToBohr);
                float y = -Ny/2-cy;
                for(size_t j=0;j<m_columns;++j)
                {
                    *y4 = _mm_set_ps1((y)*PC::AngstromToBohr);
                    *y4 = _mm_mul_ps(*z4,*y4);
                    for(size_t i=0;i<m_rows/4;++i)
                    {
                        ptr[p] = _mm_mul_ps(*y4,exp4[p]);
                        ++p;
                    }
                    y+=step;
                }
                z+=step;
            }

        #endif
    }
    else
    {
        for(size_t i=0;i<m_rows;++i)
        {
            float y = -Ny/2-cy;
            for(size_t j=0;j<m_columns;++j)
            {
                float z = -Nz/2-cz;
                for(size_t k=0;k<m_files;++k)
                {
                    (*this)(i,j,k) = exponential(i,j,k)*(y)*PC::AngstromToBohr*(z)*PC::AngstromToBohr;
                    z += step;
                }
                y += step;
            }
        }
    }
#else
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*(y)*PC::AngstromToBohr*(z)*PC::AngstromToBohr;
                z += step;
            }
            y += step;
        }
    }
#endif
}


void OrbitalArray::CalculateOrbitalFxxx( const OrbitalArray &exponential, const float cx, const float Nx, const float step)
{
#ifdef WITH_SIMD
    if  ((m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) )
    {
        __m128* ptr = (__m128*) (m_data);
        __m128* exp4 = (__m128*) (exponential.m_data);
        __m128* x4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        #ifdef ROW_MAJOR
            size_t p=0;
            float x = -Nx/2-cx;
            for(size_t i=0;i<m_rows;++i)
            {
                *x4 = _mm_set_ps1((x)*PC::AngstromToBohr);
                *x4 = _mm_mul_ps(*x4,*x4);
                *x4 = _mm_mul_ps(*x4,*x4);
                for(size_t j=0;j<m_columns;++j)
                {
                    for(size_t k=0;k<m_files/4;++k)
                    {
                        ptr[p] = _mm_mul_ps(*x4,exp4[p]);
                        ++p;
                    }
                }
                x+=step;
            }
        #else
            size_t p=0;
            for(size_t k=0;k<m_files;++k)
            {
                for(size_t j=0;j<m_columns;++j)
                {
                    float x = -Nx/2-cx;
                    for(size_t i=0;i<m_rows/4;++i)
                    {
                        *x4 = _mm_set_ps((x+3*step)*PC::AngstromToBohr,(x+2*step)*PC::AngstromToBohr,(x+step)*PC::AngstromToBohr,(x)*PC::AngstromToBohr);
                        *x4 = _mm_mul_ps(*x4,*x4);
                        *x4 = _mm_mul_ps(*x4,*x4);
                        ptr[p] = _mm_mul_ps(*x4,exp4[p]);
                        ++p;
                        x+=4*step;
                    }
                }
            }
        #endif
    }
    else
    {
        float x = -Nx/2-cx;
        for(size_t i=0;i<m_rows;++i)
        {
            for(size_t j=0;j<m_columns;++j)
            {
                for(size_t k=0;k<m_files;++k)
                {
                    (*this)(i,j,k) = exponential(i,j,k)*(x)*PC::AngstromToBohr*(x)*PC::AngstromToBohr*(x)*PC::AngstromToBohr;
                }
            }
            x += step;
        }
    }
#else
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        for(size_t j=0;j<m_columns;++j)
        {
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*(x)*PC::AngstromToBohr*(x)*PC::AngstromToBohr*(x)*PC::AngstromToBohr;
            }
        }
        x += step;
    }
#endif
}


void OrbitalArray::CalculateOrbitalFyyy( const OrbitalArray &exponential, const float cy, const float Ny, const float step)
{
#ifdef WITH_SIMD
    if  ((m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) )
    {
        __m128* ptr = (__m128*) (m_data);
        __m128* exp4 = (__m128*) (exponential.m_data);
        __m128* y4 = (__m128*) (_mm_malloc(4*sizeof(float),16));

        #ifdef ROW_MAJOR
            size_t p=0;
            for(size_t i=0;i<m_rows;++i)
            {
                float y = -Ny/2-cy;
                for(size_t j=0;j<m_columns;++j)
                {
                    *y4 = _mm_set_ps1((y)*PC::AngstromToBohr);
                    *y4 = _mm_mul_ps(*y4,*y4);
                    *y4 = _mm_mul_ps(*y4,*y4);
                    for(size_t k=0;k<m_files/4;++k)
                    {
                        ptr[p] = _mm_mul_ps(*y4,exp4[p]);
                        ++p;
                    }
                    y+=step;
                }
            }
        #else
            size_t p=0;
            for(size_t k=0;k<m_files;++k)
            {
                float y = -Ny/2-cy;
                for(size_t j=0;j<m_columns;++j)
                {
                    *y4 = _mm_set_ps1((y)*PC::AngstromToBohr);
                    *y4 = _mm_mul_ps(*y4,*y4);
                    *y4 = _mm_mul_ps(*y4,*y4);
                    for(size_t i=0;i<m_rows/4;++i)
                    {
                        ptr[p] = _mm_mul_ps(*y4,exp4[p]);
                        ++p;
                    }
                    y+=step;
                }
            }
        #endif
    }
    else
    {
        for(size_t i=0;i<m_rows;++i)
        {
            float y = -Ny/2-cy;
            for(size_t j=0;j<m_columns;++j)
            {
                for(size_t k=0;k<m_files;++k)
                {
                    (*this)(i,j,k) = exponential(i,j,k)*(y)*PC::AngstromToBohr*(y)*PC::AngstromToBohr*(y)*PC::AngstromToBohr;
                }
                y += step;
            }
        }
    }
#else
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*(y)*PC::AngstromToBohr*(y)*PC::AngstromToBohr*(y)*PC::AngstromToBohr;
            }
            y += step;
        }
    }
#endif
}


void OrbitalArray::CalculateOrbitalFzzz( const OrbitalArray &exponential, const float cz, const float Nz, const float step)
{
#ifdef WITH_SIMD
    if  ((m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) )
    {
        __m128* ptr = (__m128*) (m_data);
        __m128* exp4 = (__m128*) (exponential.m_data);
        __m128* z4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        #ifdef ROW_MAJOR
            size_t p=0;
            for(size_t i=0;i<m_rows;++i)
            {
                for(size_t j=0;j<m_columns;++j)
                {
                    float z = -Nz/2-cz;
                    for(size_t k=0;k<m_files/4;++k)
                    {
                        *z4 = _mm_set_ps((z+3*step)*PC::AngstromToBohr,(z+2*step)*PC::AngstromToBohr,(z+step)*PC::AngstromToBohr,(z)*PC::AngstromToBohr);
                        *z4 = _mm_mul_ps(*z4,*z4);
                        *z4 = _mm_mul_ps(*z4,*z4);
                        ptr[p] = _mm_mul_ps(*z4,exp4[p]);
                        ++p;
                        z+=4*step;
                    }
                }
            }
        #else
            size_t p=0;
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                *z4 = _mm_set_ps1((z)*PC::AngstromToBohr);
                *z4 = _mm_mul_ps(*z4,*z4);
                *z4 = _mm_mul_ps(*z4,*z4);
                for(size_t j=0;j<m_columns;++j)
                {
                    for(size_t i=0;i<m_rows/4;++i)
                    {
                        ptr[p] = _mm_mul_ps(*z4,exp4[p]);
                        ++p;
                    }
                }
                z+=step;
            }
        #endif
    }
    else
    {
        for(size_t i=0;i<m_rows;++i)
        {
            for(size_t j=0;j<m_columns;++j)
            {
                float z = -Nz/2-cz;
                for(size_t k=0;k<m_files;++k)
                {
                    (*this)(i,j,k) = exponential(i,j,k)*(z)*PC::AngstromToBohr*(z)*PC::AngstromToBohr*(z)*PC::AngstromToBohr;
                    z += step;
                }
            }
        }
    }
#else
    for(size_t i=0;i<m_rows;++i)
    {
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*(z)*PC::AngstromToBohr*(z)*PC::AngstromToBohr*(z)*PC::AngstromToBohr;
                z += step;
            }
        }
    }
#endif
}


void OrbitalArray::CalculateOrbitalFxxy( const OrbitalArray &exponential, const float cx, const float cy, const float Nx, const float Ny, const float step)
{
#ifdef WITH_SIMD
    if  ((m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) )
    {
        __m128* ptr = (__m128*) (m_data);
        __m128* exp4 = (__m128*) (exponential.m_data);
        __m128* x4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        __m128* y4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        #ifdef ROW_MAJOR
            size_t p=0;
            float x = -Nx/2-cx;
            for(size_t i=0;i<m_rows;++i)
            {
                *x4 = _mm_set_ps1((x)*PC::AngstromToBohr);
                *x4 = _mm_mul_ps(*x4,*x4);
                float y = -Ny/2-cy;
                for(size_t j=0;j<m_columns;++j)
                {
                    *y4 = _mm_set_ps1((y)*PC::AngstromToBohr);
                    *y4 = _mm_mul_ps(*x4,*y4);
                    for(size_t k=0;k<m_files/4;++k)
                    {
                        ptr[p] = _mm_mul_ps(*y4,exp4[p]);
                        ++p;
                    }
                    y+=step;
                }
                x+=step;
            }
        #else
            size_t p=0;
            for(size_t k=0;k<m_files;++k)
            {
                float y = -Ny/2-cy;
                for(size_t j=0;j<m_columns;++j)
                {
                    *y4 = _mm_set_ps1((y)*PC::AngstromToBohr);
                    float x = -Nx/2-cx;
                    for(size_t i=0;i<m_rows/4;++i)
                    {
                        *x4 = _mm_set_ps((x+3*step)*PC::AngstromToBohr,(x+2*step)*PC::AngstromToBohr,(x+step)*PC::AngstromToBohr,(x)*PC::AngstromToBohr);
                        *x4 = _mm_mul_ps(*x4,*x4);
                        *x4 = _mm_mul_ps(*x4,*y4);
                        ptr[p] = _mm_mul_ps(*x4,exp4[p]);
                        ++p;
                        x+=4*step;
                    }
                    y+=step;
                }
            }
        #endif
    }
    else
    {
        float x = -Nx/2-cx;
        for(size_t i=0;i<m_rows;++i)
        {
            float y = -Ny/2-cy;
            for(size_t j=0;j<m_columns;++j)
            {
                for(size_t k=0;k<m_files;++k)
                {
                    (*this)(i,j,k) = exponential(i,j,k)*(x)*PC::AngstromToBohr*(y)*PC::AngstromToBohr*(x)*PC::AngstromToBohr;
                }
                y += step;
            }
            x += step;
        }
    }
#else
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*(x)*PC::AngstromToBohr*(y)*PC::AngstromToBohr*(x)*PC::AngstromToBohr;
            }
            y += step;
        }
        x += step;
    }
#endif
}


void OrbitalArray::CalculateOrbitalFxyy( const OrbitalArray &exponential, const float cx, const float cy, const float Nx, const float Ny, const float step)
{
#ifdef WITH_SIMD
    if  ((m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) )
    {
        __m128* ptr = (__m128*) (m_data);
        __m128* exp4 = (__m128*) (exponential.m_data);
        __m128* x4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        __m128* y4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        #ifdef ROW_MAJOR
            size_t p=0;
            float x = -Nx/2-cx;
            for(size_t i=0;i<m_rows;++i)
            {
                *x4 = _mm_set_ps1((x)*PC::AngstromToBohr);
                float y = -Ny/2-cy;
                for(size_t j=0;j<m_columns;++j)
                {
                    *y4 = _mm_set_ps1((y)*PC::AngstromToBohr);
                    *y4 = _mm_mul_ps(*y4,*y4);
                    *y4 = _mm_mul_ps(*x4,*y4);
                    for(size_t k=0;k<m_files/4;++k)
                    {
                        ptr[p] = _mm_mul_ps(*y4,exp4[p]);
                        ++p;
                    }
                    y+=step;
                }
                x+=step;
            }
        #else
            size_t p=0;
            for(size_t k=0;k<m_files;++k)
            {
                float y = -Ny/2-cy;
                for(size_t j=0;j<m_columns;++j)
                {
                    *y4 = _mm_set_ps1((y)*PC::AngstromToBohr);
                    *y4 = _mm_mul_ps(*y4,*y4);
                    float x = -Nx/2-cx;
                    for(size_t i=0;i<m_rows/4;++i)
                    {
                        *x4 = _mm_set_ps((x+3*step)*PC::AngstromToBohr,(x+2*step)*PC::AngstromToBohr,(x+step)*PC::AngstromToBohr,(x)*PC::AngstromToBohr);
                        *x4 = _mm_mul_ps(*x4,*y4);
                        ptr[p] = _mm_mul_ps(*x4,exp4[p]);
                        ++p;
                        x+=4*step;
                    }
                    y+=step;
                }
            }
        #endif
    }
    else
    {
        float x = -Nx/2-cx;
        for(size_t i=0;i<m_rows;++i)
        {
            float y = -Ny/2-cy;
            for(size_t j=0;j<m_columns;++j)
            {
                for(size_t k=0;k<m_files;++k)
                {
                    (*this)(i,j,k) = exponential(i,j,k)*(x)*PC::AngstromToBohr*(y)*PC::AngstromToBohr*(y)*PC::AngstromToBohr;
                }
                y += step;
            }
            x += step;
        }
    }
#else
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*(x)*PC::AngstromToBohr*(y)*PC::AngstromToBohr*(y)*PC::AngstromToBohr;
            }
            y += step;
        }
        x += step;
    }
#endif
}


void OrbitalArray::CalculateOrbitalFxxz( const OrbitalArray &exponential, const float cx, const float cz, const float Nx, const float Nz, const float step)
{
#ifdef WITH_SIMD
    if  ((m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) )
    {
        __m128* ptr = (__m128*) (m_data);
        __m128* exp4 = (__m128*) (exponential.m_data);
        __m128* x4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        __m128* z4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        #ifdef ROW_MAJOR
            size_t p=0;
            float x = -Nx/2-cx;
            for(size_t i=0;i<m_rows;++i)
            {
                *x4 = _mm_set_ps1((x)*PC::AngstromToBohr);
                *x4 = _mm_mul_ps(*x4,*x4);
                for(size_t j=0;j<m_columns;++j)
                {
                    float z = -Nz/2-cz;
                    for(size_t k=0;k<m_files/4;++k)
                    {
                        *z4 = _mm_set_ps((z+3*step)*PC::AngstromToBohr,(z+2*step)*PC::AngstromToBohr,(z+step)*PC::AngstromToBohr,(z)*PC::AngstromToBohr);
                        *z4 = _mm_mul_ps(*x4,*z4);
                        ptr[p] = _mm_mul_ps(*z4,exp4[p]);
                        ++p;
                        z+=4*step;
                    }
                }
                x+=step;
            }
        #else
            size_t p=0;
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                *z4 = _mm_set_ps1((z)*PC::AngstromToBohr);
                for(size_t j=0;j<m_columns;++j)
                {
                    float x = -Nx/2-cx;
                    for(size_t i=0;i<m_rows/4;++i)
                    {
                        *x4 = _mm_set_ps((x+3*step)*PC::AngstromToBohr,(x+2*step)*PC::AngstromToBohr,(x+step)*PC::AngstromToBohr,(x)*PC::AngstromToBohr);
                        *x4 = _mm_mul_ps(*x4,*x4);
                        *x4 = _mm_mul_ps(*x4,*z4);
                        ptr[p] = _mm_mul_ps(*x4,exp4[p]);
                        ++p;
                        x+=4*step;
                    }
                }
                z+=step;
            }
        #endif
    }
    else
    {
        float x = -Nx/2-cx;
        for(size_t i=0;i<m_rows;++i)
        {
            for(size_t j=0;j<m_columns;++j)
            {
                float z = -Nz/2-cz;
                for(size_t k=0;k<m_files;++k)
                {
                    (*this)(i,j,k) = exponential(i,j,k)*(x)*PC::AngstromToBohr*(x)*PC::AngstromToBohr*(z)*PC::AngstromToBohr;
                    z += step;
                }
            }
            x += step;
        }
    }
#else
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*(x)*PC::AngstromToBohr*(x)*PC::AngstromToBohr*(z)*PC::AngstromToBohr;
                z += step;
            }
        }
        x += step;
    }
#endif
}


void OrbitalArray::CalculateOrbitalFxzz( const OrbitalArray &exponential, const float cx, const float cz, const float Nx, const float Nz, const float step)
{
#ifdef WITH_SIMD
    if  ((m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) )
    {
        __m128* ptr = (__m128*) (m_data);
        __m128* exp4 = (__m128*) (exponential.m_data);
        __m128* x4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        __m128* z4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        #ifdef ROW_MAJOR
            size_t p=0;
            float x = -Nx/2-cx;
            for(size_t i=0;i<m_rows;++i)
            {
                *x4 = _mm_set_ps1((x)*PC::AngstromToBohr);
                for(size_t j=0;j<m_columns;++j)
                {
                    float z = -Nz/2-cz;
                    for(size_t k=0;k<m_files/4;++k)
                    {
                        *z4 = _mm_set_ps((z+3*step)*PC::AngstromToBohr,(z+2*step)*PC::AngstromToBohr,(z+step)*PC::AngstromToBohr,(z)*PC::AngstromToBohr);
                        *z4 = _mm_mul_ps(*z4,*z4);
                        *z4 = _mm_mul_ps(*x4,*z4);
                        ptr[p] = _mm_mul_ps(*z4,exp4[p]);
                        ++p;
                        z+=4*step;
                    }
                }
                x+=step;
            }
        #else
            size_t p=0;
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                *z4 = _mm_set_ps1((z)*PC::AngstromToBohr);
                *z4 = _mm_mul_ps(*z4,*z4);
                for(size_t j=0;j<m_columns;++j)
                {
                    float x = -Nx/2-cx;
                    for(size_t i=0;i<m_rows/4;++i)
                    {
                        *x4 = _mm_set_ps((x+3*step)*PC::AngstromToBohr,(x+2*step)*PC::AngstromToBohr,(x+step)*PC::AngstromToBohr,(x)*PC::AngstromToBohr);
                        *x4 = _mm_mul_ps(*x4,*z4);
                        ptr[p] = _mm_mul_ps(*x4,exp4[p]);
                        ++p;
                        x+=4*step;
                    }
                }
                z+=step;
            }
        #endif
    }
    else
    {
        float x = -Nx/2-cx;
        for(size_t i=0;i<m_rows;++i)
        {
            for(size_t j=0;j<m_columns;++j)
            {
                float z = -Nz/2-cz;
                for(size_t k=0;k<m_files;++k)
                {
                    (*this)(i,j,k) = exponential(i,j,k)*(x)*PC::AngstromToBohr*(z)*PC::AngstromToBohr*(z)*PC::AngstromToBohr;
                    z += step;
                }
            }
            x += step;
        }
    }
#else
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*(x)*PC::AngstromToBohr*(z)*PC::AngstromToBohr*(z)*PC::AngstromToBohr;
                z += step;
            }
        }
        x += step;
    }
#endif
}


void OrbitalArray::CalculateOrbitalFyyz( const OrbitalArray &exponential, const float cy, const float cz, const float Ny, const float Nz, const float step)
{
#ifdef WITH_SIMD
    if  ((m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) )
    {
        __m128* ptr = (__m128*) (m_data);
        __m128* exp4 = (__m128*) (exponential.m_data);
        __m128* y4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        __m128* z4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        #ifdef ROW_MAJOR
            size_t p=0;
            for(size_t i=0;i<m_rows;++i)
            {
                float y = -Ny/2-cy;
                for(size_t j=0;j<m_columns;++j)
                {
                    *y4 = _mm_set_ps1((y)*PC::AngstromToBohr);
                    *y4 = _mm_mul_ps(*y4,*y4);
                    float z = -Nz/2-cz;
                    for(size_t k=0;k<m_files/4;++k)
                    {
                        *z4 = _mm_set_ps((z+3*step)*PC::AngstromToBohr,(z+2*step)*PC::AngstromToBohr,(z+step)*PC::AngstromToBohr,(z)*PC::AngstromToBohr);
                        *z4 = _mm_mul_ps(*z4,*y4);
                        ptr[p] = _mm_mul_ps(*z4,exp4[p]);
                        ++p;
                        z+=4*step;
                    }
                    y+=step;
                }
            }
        #else
            size_t p=0;
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                *z4 = _mm_set_ps1((z)*PC::AngstromToBohr);
                float y = -Ny/2-cy;
                for(size_t j=0;j<m_columns;++j)
                {
                    *y4 = _mm_set_ps1((y)*PC::AngstromToBohr);
                    *y4 = _mm_mul_ps(*y4,*y4);
                    *y4 = _mm_mul_ps(*z4,*y4);
                    for(size_t i=0;i<m_rows/4;++i)
                    {
                        ptr[p] = _mm_mul_ps(*y4,exp4[p]);
                        ++p;
                    }
                    y+=step;
                }
                z+=step;
            }

        #endif
    }
    else
    {
        for(size_t i=0;i<m_rows;++i)
        {
            float y = -Ny/2-cy;
            for(size_t j=0;j<m_columns;++j)
            {
                float z = -Nz/2-cz;
                for(size_t k=0;k<m_files;++k)
                {
                    (*this)(i,j,k) = exponential(i,j,k)*(y)*PC::AngstromToBohr*(y)*PC::AngstromToBohr*(z)*PC::AngstromToBohr;
                    z += step;
                }
                y += step;
            }
        }
    }
#else
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*(y)*PC::AngstromToBohr*(y)*PC::AngstromToBohr*(z)*PC::AngstromToBohr;
                z += step;
            }
            y += step;
        }
    }
#endif
}


void OrbitalArray::CalculateOrbitalFyzz( const OrbitalArray &exponential, const float cy, const float cz, const float Ny, const float Nz, const float step)
{
#ifdef WITH_SIMD
    if  ((m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) )
    {
        __m128* ptr = (__m128*) (m_data);
        __m128* exp4 = (__m128*) (exponential.m_data);
        __m128* y4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        __m128* z4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        #ifdef ROW_MAJOR
            size_t p=0;
            for(size_t i=0;i<m_rows;++i)
            {
                float y = -Ny/2-cy;
                for(size_t j=0;j<m_columns;++j)
                {
                    *y4 = _mm_set_ps1((y)*PC::AngstromToBohr);
                    float z = -Nz/2-cz;
                    for(size_t k=0;k<m_files/4;++k)
                    {
                        *z4 = _mm_set_ps((z+3*step)*PC::AngstromToBohr,(z+2*step)*PC::AngstromToBohr,(z+step)*PC::AngstromToBohr,(z)*PC::AngstromToBohr);
                        *z4 = _mm_mul_ps(*z4,*z4);
                        *z4 = _mm_mul_ps(*z4,*y4);
                        ptr[p] = _mm_mul_ps(*z4,exp4[p]);
                        ++p;
                        z+=4*step;
                    }
                    y+=step;
                }
            }
        #else
            size_t p=0;
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files/4;++k)
            {
                *z4 = _mm_set_ps1((z)*PC::AngstromToBohr);
                *z4 = _mm_mul_ps(*z4,*z4);
                float y = -Ny/2-cy;
                for(size_t j=0;j<m_columns;++j)
                {
                    *y4 = _mm_set_ps1((y)*PC::AngstromToBohr);
                    *y4 = _mm_mul_ps(*z4,*y4);
                    for(size_t i=0;i<m_rows;++i)
                    {
                        ptr[p] = _mm_mul_ps(*y4,exp4[p]);
                        ++p;
                    }
                    y+=step;
                }
                z+=step;
            }

        #endif
    }
    else
    {
        for(size_t i=0;i<m_rows;++i)
        {
            float y = -Ny/2-cy;
            for(size_t j=0;j<m_columns;++j)
            {
                float z = -Nz/2-cz;
                for(size_t k=0;k<m_files;++k)
                {
                    (*this)(i,j,k) = exponential(i,j,k)*(z)*PC::AngstromToBohr*(y)*PC::AngstromToBohr*(z)*PC::AngstromToBohr;
                    z += step;
                }
                y += step;
            }
        }
    }
#else
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*(z)*PC::AngstromToBohr*(y)*PC::AngstromToBohr*(z)*PC::AngstromToBohr;
                z += step;
            }
            y += step;
        }
    }
#endif
}


void OrbitalArray::CalculateOrbitalFxyz( const OrbitalArray &exponential, const float cx, const float cy, const float cz, const float Nx, const float Ny, const float Nz, const float step)
{
#ifdef WITH_SIMD
    if  ((m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) )
    {
        __m128* ptr = (__m128*) (m_data);
        __m128* exp4 = (__m128*) (exponential.m_data);
        __m128* x4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        __m128* y4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        __m128* z4 = (__m128*) (_mm_malloc(4*sizeof(float),16));
        #ifdef ROW_MAJOR
            size_t p=0;
            float x = -Nx/2-cx;
            for(size_t i=0;i<m_rows;++i)
            {
                *x4 = _mm_set_ps1((x)*PC::AngstromToBohr);
                float y = -Ny/2-cy;
                for(size_t j=0;j<m_columns;++j)
                {
                    *y4 = _mm_set_ps1((y)*PC::AngstromToBohr);
                    *y4 = _mm_mul_ps(*x4,*y4);
                    float z = -Nz/2-cz;
                    for(size_t k=0;k<m_files/4;++k)
                    {
                        *z4 = _mm_set_ps((z+3*step)*PC::AngstromToBohr,(z+2*step)*PC::AngstromToBohr,(z+step)*PC::AngstromToBohr,(z)*PC::AngstromToBohr);
                        *z4 = _mm_mul_ps(*y4,*z4);
                        ptr[p] = _mm_mul_ps(*z4,exp4[p]);
                        ++p;
                        z+=4*step;
                    }
                    y+=step;
                }
                x+=step;
            }
        #else
            size_t p=0;
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                *z4 = _mm_set_ps1((z)*PC::AngstromToBohr);
                float y = -Ny/2-cy;
                for(size_t j=0;j<m_columns;++j)
                {
                    *y4 = _mm_set_ps1((y)*PC::AngstromToBohr);
                    *y4 = _mm_mul_ps(*z4,*y4);
                    float x = -Nx/2-cx;
                    for(size_t i=0;i<m_rows/4;++i)
                    {
                        *x4 = _mm_set_ps((x+3*step)*PC::AngstromToBohr,(x+2*step)*PC::AngstromToBohr,(x+step)*PC::AngstromToBohr,(x)*PC::AngstromToBohr);
                        *x4 = _mm_mul_ps(*x4,*y4);
                        ptr[p] = _mm_mul_ps(*x4,exp4[p]);
                        ++p;
                        x+=4*step;
                    }
                    y+=step;
                }
                z+=step;
            }
        #endif
    }
    else
    {
        float x = -Nx/2-cx;
        for(size_t i=0;i<m_rows;++i)
        {
            float y = -Ny/2-cy;
            for(size_t j=0;j<m_columns;++j)
            {
                float z = -Nz/2-cz;
                for(size_t k=0;k<m_files;++k)
                {
                    (*this)(i,j,k) = exponential(i,j,k)*(x)*PC::AngstromToBohr*(y)*PC::AngstromToBohr*(z)*PC::AngstromToBohr;
                    z += step;
                }
                y += step;
            }
            x += step;
        }
    }
#else
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*(x)*PC::AngstromToBohr*(y)*PC::AngstromToBohr*(z)*PC::AngstromToBohr;
                z += step;
            }
            y += step;
        }
        x += step;
    }
#endif
}


void OrbitalArray::CalculateOrbitalDY0( const OrbitalArray &exponential, const float cx, const float cy, const float cz, const float Nx, const float Ny, const float Nz, const float step)
{
#ifdef WITH_SIMD
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstDY0*(-(x)*(x)-(y)*(y)+2*(z)*(z))*PC::AngstromToBohr*PC::AngstromToBohr;
                z += step;
            }
            y += step;
        }
        x += step;
    }
#else
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstDY0*(-(x)*(x)-(y)*(y)+2*(z)*(z))*PC::AngstromToBohr*PC::AngstromToBohr;
                z += step;
            }
            y += step;
        }
        x += step;
    }
#endif
}


void OrbitalArray::CalculateOrbitalDY1( const OrbitalArray &exponential, const float cx, const float cz, const float Nx, const float Nz, const float step)
{
#ifdef WITH_SIMD
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstDY1*(x)*(z)*PC::AngstromToBohr*PC::AngstromToBohr;
                z += step;
            }
        }
        x += step;
    }
#else
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstDY1*(x)*(z)*PC::AngstromToBohr*PC::AngstromToBohr;
                z += step;
            }
        }
        x += step;
    }
#endif
}


void OrbitalArray::CalculateOrbitalDY2( const OrbitalArray &exponential, const float cy, const float cz, const float Ny, const float Nz, const float step)
{
#ifdef WITH_SIMD
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstDY2*(y)*(z)*PC::AngstromToBohr*PC::AngstromToBohr;
                z += step;
            }
            y += step;
        }
    }
#else
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstDY2*(y)*(z)*PC::AngstromToBohr*PC::AngstromToBohr;
                z += step;
            }
            y += step;
        }
    }
#endif
}


void OrbitalArray::CalculateOrbitalDY3( const OrbitalArray &exponential, const float cx, const float cy, const float Nx, const float Ny, const float step)
{
#ifdef WITH_SIMD
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstDY3*((x)*(x)-(y)*(y))*PC::AngstromToBohr*PC::AngstromToBohr;
            }
            y += step;
        }
        x += step;
    }
#else
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstDY3*((x)*(x)-(y)*(y))*PC::AngstromToBohr*PC::AngstromToBohr;
            }
            y += step;
        }
        x += step;
    }
#endif
}


void OrbitalArray::CalculateOrbitalDY4( const OrbitalArray &exponential, const float cx, const float cy, const float Nx, const float Ny, const float step)
{
#ifdef WITH_SIMD
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstDY4*((x)*(y))*PC::AngstromToBohr*PC::AngstromToBohr;
            }
            y += step;
        }
        x += step;
    }
#else
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstDY4*((x)*(y))*PC::AngstromToBohr*PC::AngstromToBohr;
            }
            y += step;
        }
        x += step;
    }
#endif
}


void OrbitalArray::CalculateOrbitalFY0( const OrbitalArray &exponential, const float cx, const float cy, const float cz, const float Nx, const float Ny, const float Nz, const float step)
{
#ifdef WITH_SIMD
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstFY0*(2*(z)*(z)-3*(x)*(x)-3*(y)*(y))*(z)*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
                z += step;
            }
            y += step;
        }
        x += step;
    }
#else
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstFY0*(2*(z)*(z)-3*(x)*(x)-3*(y)*(y))*(z)*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
                z += step;
            }
            y += step;
        }
        x += step;
    }
#endif
}


void OrbitalArray::CalculateOrbitalFY1( const OrbitalArray &exponential, const float cx, const float cy, const float cz, const float Nx, const float Ny, const float Nz, const float step)
{
#ifdef WITH_SIMD
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstFY1*(4*(z)*(z)-(x)*(x)-(y)*(y))*(x)*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
                z += step;
            }
            y += step;
        }
        x += step;
    }
#else
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstFY1*(4*(z)*(z)-(x)*(x)-(y)*(y))*(x)*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
                z += step;
            }
            y += step;
        }
        x += step;
    }
#endif
}


void OrbitalArray::CalculateOrbitalFY2( const OrbitalArray &exponential, const float cx, const float cy, const float cz, const float Nx, const float Ny, const float Nz, const float step)
{
#ifdef WITH_SIMD
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstFY2*(4*(z)*(z)-(x)*(x)-(y)*(y))*(y)*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
                z += step;
            }
            y += step;
        }
        x += step;
    }
#else
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstFY2*(4*(z)*(z)-(x)*(x)-(y)*(y))*(y)*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
                z += step;
            }
            y += step;
        }
        x += step;
    }
#endif
}


void OrbitalArray::CalculateOrbitalFY3( const OrbitalArray &exponential, const float cx, const float cy, const float cz, const float Nx, const float Ny, const float Nz, const float step)
{
#ifdef WITH_SIMD
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstFY3*((x)*(x)-(y)*(y))*(z)*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
                z += step;
            }
            y += step;
        }
        x += step;
    }
#else
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstFY3*((x)*(x)-(y)*(y))*(z)*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
                z += step;
            }
            y += step;
        }
        x += step;
    }
#endif
}


void OrbitalArray::CalculateOrbitalFY4( const OrbitalArray &exponential, const float cx, const float cy, const float cz, const float Nx, const float Ny, const float Nz, const float step)
{
#ifdef WITH_SIMD
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstFY4*(x)*(y)*(z)*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
                z += step;
            }
            y += step;
        }
        x += step;
    }
#else
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            float z = -Nz/2-cz;
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstFY4*(x)*(y)*(z)*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
                z += step;
            }
            y += step;
        }
        x += step;
    }
#endif
}


void OrbitalArray::CalculateOrbitalFY5( const OrbitalArray &exponential, const float cx, const float cy, const float Nx, const float Ny, const float step)
{
#ifdef WITH_SIMD
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstFY5*((x)*(x)-3*(y)*(y))*(x)*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
            }
            y += step;
        }
        x += step;
    }
#else
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstFY5*((x)*(x)-3*(y)*(y))*(x)*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
            }
            y += step;
        }
        x += step;
    }
#endif
}


void OrbitalArray::CalculateOrbitalFY6( const OrbitalArray &exponential, const float cx, const float cy, const float Nx, const float Ny, const float step)
{
#ifdef WITH_SIMD
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstFY6*(3*(x)*(x)-(y)*(y))*(y)*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
            }
            y += step;
        }
        x += step;
    }
#else
    float x = -Nx/2-cx;
    for(size_t i=0;i<m_rows;++i)
    {
        float y = -Ny/2-cy;
        for(size_t j=0;j<m_columns;++j)
        {
            for(size_t k=0;k<m_files;++k)
            {
                (*this)(i,j,k) = exponential(i,j,k)*PC::ConstFY6*(3*(x)*(x)-(y)*(y))*(y)*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
            }
            y += step;
        }
        x += step;
    }
#endif
}
