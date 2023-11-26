/*****************************************************************************************
                            fidarray.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <iostream>
#include "fidarray.h"

fidarray::fidarray()
{
}

fidarray::fidarray(size_t n) : std::valarray< std::complex<float> >(n)
{}
fidarray::~fidarray()
{}

fidarray& fidarray::shift(int np)
{
    if( np == 0 ) return *this;
    fidarray tmp(this->size());
    tmp=(*this);
    if( np < 0)
    {
        np*=-1;
        size_t i,j;
        for(i=0;i< (size_t) np;i++)
            (*this)[i]=0;
        for(i=(size_t) np,j=0;i<this->size();i++,j++)
            (*this)[i]=tmp[j];
    }

    if(np > 0)
    {
        size_t i,j;
        for(i=0,j=(size_t) np;i< this->size()-(size_t) np;i++,j++)
        {
            (*this)[i]=tmp[j];
        }
        for(i=this->size()-(size_t)np;i< this->size();i++)
        {
            (*this)[i]=0;
        }
    }

    return *this;

}

fidarray& fidarray::cshift(int np)
{
    if( np == 0 ) return *this;
    fidarray tmp(this->size());
    tmp=(*this);
    if( np < 0)
    {
        np*=-1;
        size_t i,j;
        for(i=0,j=(this->size()-(size_t)np);i< (size_t)np;i++,j++)
            (*this)[i]=tmp[j];
        for(i=(size_t)np,j=0;i<this->size();i++,j++)
            (*this)[i]=tmp[j];
    }

    if(np > 0)
    {
        size_t i,j;
        for(i=0,j=(size_t)np;i< this->size()-(size_t)np;i++,j++)
        {
            (*this)[i]=tmp[j];
        }
        for(i=this->size()-(size_t)np,j=0;i< this->size();i++,j++)
        {
            (*this)[i]=tmp[j];
        }
    }

    return *this;
}

void fidarray::grow(size_t n)
{

    if ( n==size() ) return;

    if( n>size() )
    {
        fidarray tmp(size());
        std::slice s(0,size(),1);
        std::slice s1(size(),n-size(),1);
        tmp=*this;
        resize(n,std::complex<float>(0,0));
    //    (*this)[s]=tmp;
       (*this)[s]=tmp[s];
       // k=tmp[s];
       // const_cast< std::slice_array< std::complex<float> >&  >((*this)[s]);
        (*this)[s1]=std::complex<float>(0,0);

    }
    else
    {
        fidarray tmp(n);
        tmp=(*this);
        resize(n);
        (*this)=tmp;
    }
}

fidarray& fidarray::quadrature()
{

    for (size_t i=1; i< size(); i+=2)
    {
        (*this)[i]*=-1.0f;
    }

    return *this;


}

#include "gsl_wavelet.h"
void fidarray::wlt(int sense, int coefficients)
{
    if( coefficients != 4 &&  coefficients != 12 && coefficients != 20)
    {
        std::cerr << "No data for Daubechies(" << coefficients << ")" << std::endl;
        return;
    }
    size_t n=this->size();
    gsl_wavelet *w = gsl_wavelet_alloc (gsl_wavelet_daubechies, coefficients);
    gsl_wavelet_workspace *work = gsl_wavelet_workspace_alloc (n);

    float* buffer=(float*)(&((*this)[0]));
    if (sense == 1 )
    {
        gsl_wavelet_transform_forward(w, buffer , 2, n, work);
        gsl_wavelet_transform_forward(w, buffer+1, 2, n, work);}

    else if ( sense == -1 )
    {
        gsl_wavelet_transform_inverse(w, buffer, 2, n, work); gsl_wavelet_transform_inverse(w, buffer+1, 2, n, work);
    }

}

void fidarray::fft(int sense)
{

    fidarray& input=(*this);
    int result ;
    size_t dual;
    size_t bit;
    size_t logn = 0;

    size_t n=input.size();

    if (n == 1) /* identity operation */
    {
        return  ;
    }

    /* make sure that n is a power of 2 */

    bool b;
    result = logsize(b);

    if (!b)
    {
        std::cerr << "fidarray::fft not a power of two" <<std::endl;
        return;
    }
    else
    {
        logn = result ;

    }

    /* bit reverse the ordering of input data for decimation in time algorithm */

    input.bitreverse();

    /* apply fft recursion */

    dual = 1;

    for (bit = 0; bit < logn; bit++)
    {
        std::complex<float> w(1.0,0.0);

        const double theta = 2.0 * ((int) sense) * M_PI / (2.0 * (double) dual);

        double t = sin (theta / 2.0);
        std::complex<float> s(-2.0*t*t,sin(theta));

        size_t a, b;

        /* a = 0 */

        for (b = 0; b < n; b += 2 * dual)
        {
            const size_t i = b ;
            const size_t j = b + dual;

            std::complex<float> z1=input[j];
            std::complex<float> wd=z1;
            input[j]=input[i]-wd;
            input[i]+=wd;
        }

        /* a = 1 .. (dual-1) */

        for (a = 1; a < dual; a++)
        {

            /* trignometric recurrence for w-> exp(i theta) w */

            {
                std::complex<float> tmp=w+s*w;
                w=tmp;

            }

            for (b = 0; b < n; b += 2 * dual)
            {
                const size_t i = b + a;
                const size_t j = b + a + dual;
                std::complex<float> z1=input[j];;
                std::complex<float> wd=w*z1;

                input[j]=input[i]-wd;
                input[i]+=wd;

            }
        }
        dual *= 2;
    }
}

void fidarray::clear()
{
    resize(0);
}

float fidarray::sum(bool real)
{
    float suma=0.0f;

    if(real)
        for(register int i=0;i<(int)size();i++)
        {
            suma+=(*this)[i].real();
        }

    else
        for(register int i=0;i<(int)size();i++)
        {
            suma+=(*this)[i].imag();
        }

    return suma;
}

void fidarray::magnitude()
{
#warning check
    for(int i=0;i<(int)size();++i)
    {
        (*this)[i]=std::norm((*this)[i]);
    }
}

void fidarray::power()
{
#warning check
    for(int i=0;i<(int)size();++i)
    {
        (*this)[i]=std::abs((*this)[i]);;
    }

}

void fidarray::bitreverse()
{
    size_t i;
    size_t j = 0;

    size_t n=size();

    for (i = 0; i < n - 1; i++)
    {
        size_t k = n / 2 ;

        if (i < j)
            std::swap((*this)[j],(*this)[i]);

        while (k <= j)
        {
            j = j - k ;
            k = k / 2 ;
        }

        j += k ;
    }

    return;
}

void fidarray::bitreverse_real()
{
    /* This is the Goldrader bit-reversal algorithm */

    size_t i;
    size_t j = 0;

    size_t n=2*(this->size());
    float* input=(float*)&((*this)[0]);

    // logn = 0 ; /* not needed for this algorithm */

    for (i = 0; i < n - 1; i++)
    {
        size_t k = n / 2 ;

        if (i < j)
        {
            std::swap(input[i],input[j]);
            /* const BASE tmp = VECTOR(data,stride,i);
       VECTOR(data,stride,i) = VECTOR(data,stride,j);
       VECTOR(data,stride,j) = tmp;  */
        }

        while (k <= j)
        {
            j = j - k ;
            k = k / 2 ;
        }

        j += k ;
    }


}
size_t fidarray::logsize(bool& ispoweroftwo)
{
    size_t n=size();
    size_t ntest ;
    size_t binary_logn = 0 ;
    size_t k = 1;

    while (k < n)
    {
        k *= 2;
        binary_logn++;
    }

    ntest = (1 << binary_logn) ;

    if (n != ntest )
    {
        ispoweroftwo=false;
    }
    else ispoweroftwo=true;

    return binary_logn;
}



void fidarray::fft_real()
{

    int result ;
    size_t p, p_1, q;
    size_t i;
    size_t logn = 0;
    size_t n=2*this->size();
    if (n == 1) /* identity operation */
    {
        return;
    }

    /* make sure that n is a power of 2 */
    bool b;
    result = logsize(b) ;

    if (!b)
    {
        std::cerr << "fidarray::fft_real n is not a power of 2"<<std::endl;
        return;
    }
    else
    {
        logn = 2*result ;
    }

    /* bit reverse the ordering of input data for decimation in time algorithm */

    this->bitreverse_real();
    float* input=(float*)&((*this)[0]);
    // apply fft recursion

    p = 1; q = n ;

    // timeval start,end;

    for (i = 1; i <= logn; i++)
    {
        size_t a, b;

        p_1 = p ;
        p = 2 * p ;
        q = q / 2 ;

        // a = 0

        for (b = 0; b < q; b++)
        {
            double t0_real = input[b*p]+input[b*p+p_1];
            double t1_real = input[b*p]-input[b*p+p_1];

            input[b*p]=t0_real;
            input[b*p+p_1]=t1_real;
        }

        // a = 1 ... p_{i-1}/2 - 1

        {
            double w_real = 1.0;
            double w_imag = 0.0;

            const double theta = - 2.0 * M_PI / p;

            const double s = sin (theta);
            const double t = sin (theta / 2.0);
            const double s2 = 2.0 * t * t;

            for (a = 1; a < (p_1)/2; a++)
            {
                // trignometric recurrence for w-> exp(i theta) w

                {
                    const double tmp_real = w_real - s * w_imag - s2 * w_real;
                    const double tmp_imag = w_imag + s * w_real - s2 * w_imag;
                    w_real = tmp_real;
                    w_imag = tmp_imag;
                }

                size_t counter=0;



                for (b = 0; b < q; b++)
                {
                    int c[4];
                    c[0]=counter+a;
                    c[1]=counter+p_1-a;
                    c[2]=counter+p_1+a;
                    c[3]=counter+p-a;
                    double z0_real = input[c[0]] ;
                    double z0_imag = input[c[1]] ;
                    double z1_real = input[c[2]] ;
                    double z1_imag = input[c[3]] ;

                    // t0 = z0 + w * z1

                    double t0_real = z0_real + w_real * z1_real - w_imag * z1_imag;
                    double t0_imag = z0_imag + w_real * z1_imag + w_imag * z1_real;

                    // t1 = z0 - w * z1

                    double t1_real = z0_real - w_real * z1_real + w_imag * z1_imag;
                    double t1_imag = z0_imag - w_real * z1_imag - w_imag * z1_real;

                    input[ c[0]] = t0_real ;
                    input[c[3]] = t0_imag ;

                    input[c[1]] = t1_real ;
                    input[c[2]] = -t1_imag ;

                    counter+=p;
                }
            }


        }

        if (p_1 >  1)
        {
            for (b = 0; b < q; b++)
            {
                // a = p_{i-1}/2

                input[b*p + p - p_1/2] *= -1 ;
            }
        }
    }

    return;
}

void fidarray::real2complex()
{

    fidarray temp;
    size_t n=this->size();
    temp.grow(n/2);
    size_t i;
    for(i=0;i<n/2;i++)
    {
        temp[i].real((*this)[i].real());
        temp[i].imag((*this)[n-i].imag());
    }

    this->grow(n/2);
    (*this)=temp;


}

void fidarray::normalize()
{
    size_t n=this->size();
    for(size_t i=0;i<n;i++)
    {
        (*this)[i]/=n;
    }
}

