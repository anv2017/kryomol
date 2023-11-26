/*****************************************************************************************
                            mathtools.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef MATHTOOLS_H
#define MATHTOOLS_H

#include <cstring>
#include <assert.h>
#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <typeinfo>
#include "exception.h"
#include "toolsexport.h"
#include "arrays_defines.h"
#include "physicalconstants.h"

#include <QDebug>

#ifdef WITH_MKL
#include <mkl.h>
#ifdef __GNUC__
#warning compiling with mkl, column major storage
#endif
#endif

#ifdef WITH_LAPACK
#warning compiling with lapack, column major storage
#endif

#ifdef WITH_VECLIB
#include <Accelerate/Accelerate.h>
#warning compiling with vecLib, column major storage
#endif

#ifdef WITH_CUDA
#include "cudamathtools.h"
#endif

#ifdef WITH_SIMD
#include "sse_mathfun.h"
#endif


/** */
namespace kryomol
{

/** @brief monodimensional array

     abstraction of memory-contiguos monodimensional arrays
     */

template <class X>
class D1Array
{

private:
    X* m_data;
public:
    /** Build empty array*/
    D1Array() : m_size ( 0 )
    {
        m_data=NULL;
    }
    /** Build an array of size i and undefined values*/
    explicit D1Array ( size_t i ) : m_size ( i )
    {
        m_data= new X[i];
    }
    /** Build an array of size i, with an x value for all elements*/
    D1Array ( size_t i,const X& x );
    /** copy constructor*/
    D1Array ( const D1Array& x );
    /** resize an existing array to size i, all values are undefined and data are lost*/
    void Initialize ( size_t i )
    {
        if ( m_data )
            delete [] m_data;
        m_data = new X[i];
        m_size=i;
    }
    /** resize an existing array to size i, all values are set to x and previous data are deleted*/
    void Initialize ( size_t i, const X& x );
    /** destructor all date are deleted from the heap*/
    ~D1Array()
    {
        if ( m_data )
            delete [] m_data;
    }
    void Clear()
    {
        if ( m_data )
            delete [] m_data;
        m_data=NULL;
        m_size=0;
    }
    /** @return a reference to element i*/
    X& operator () ( size_t i );
    /** @return a const reference to element i*/
    const X& operator() ( size_t i ) const;
    /** casting operator to X*, allow easy pass to math libraries*/
    operator X*()
    {
        return m_data;
    }
    /** casting operator to const X*, allow easy pass to math libraries*/
    operator const X*() const
    {
        return m_data;
    }
    /** asignment operator*/
    D1Array& operator= ( const D1Array& x );
    /** @return the size of the array*/
    size_t size() const
    {
        return m_size;
    }

private:
    size_t m_size;

};


template<class X>
D1Array<X>::D1Array ( size_t i, const X& x ) : m_size ( i )
{
    m_data = new X[i];
    for ( size_t k=0;k<i;k++ )
        m_data[k]=x;
}
template <class X>
D1Array<X>::D1Array ( const D1Array<X>& x )
{

    m_size=x.size();
    if ( x.m_data !=NULL )
    {
        m_data= new X[m_size];

        memcpy ( m_data,x,m_size*sizeof ( X ) );
    }
    else m_data = NULL;

}

template <class X>
D1Array<X>& D1Array<X>::operator= ( const D1Array<X>& x )
{
    if ( m_data )
        delete [] m_data;
    m_size=x.size();
    if ( x.m_data != NULL )
    {
        m_data= new X[m_size];
        memcpy ( m_data,x,m_size*sizeof ( X ) );
    }
    else m_data = NULL;
    return *this;

}
template<class X>
void D1Array<X>::Initialize ( size_t i, const X& x )
{
    if ( m_data )
        delete [] m_data;
    m_size=i;
    m_data = new X[i];
    for ( size_t k=0;k<i;k++ )
        m_data[k]=x;
}

template <class X>
X& D1Array<X>::operator () ( size_t i )
{
    assert ( i < m_size );
    return m_data[i];
}

template <class X>
const X& D1Array<X>::operator () ( size_t i ) const
{
    assert ( i < m_size );
    return m_data[i];
}


/** @brief bidimensional array

     abstraction of memory-contiguous bidimensional arrays
     */
template<class X>
class D2Array
{
private:
    X* m_data;
public:
    /** Build an empty N-dimensional array */
    D2Array()
    {
        m_rows=m_columns=0;
        m_data=NULL;
    }
    /** Build an array with  i rows and j columns*/
    explicit D2Array ( size_t i, size_t j ) : m_rows ( i ), m_columns ( j )
    {
        m_data = new X[i*j];
    }
    /** Build an array with  i rows and j columns*, with values x*/
    D2Array ( size_t i,size_t j, const X& x );
    /** copy constructor*/
    D2Array ( const D2Array& x );
    /** assignment operator*/
    D2Array& operator= ( const D2Array& x );
    /** sum operator */
    void operator+=(const D2Array& x);
    /** difference operator*/
    void operator-=(const D2Array& x);
    /** multiply by scalar x*/
    void operator*=(const X& x);
    /** resize an existing array to i rows and j columns, erase all previous data*/
    void Initialize ( size_t i, size_t j )
    {
        assert( i > 0 && j > 0 );
        m_rows=i;
        m_columns=j;
        if ( m_data )
            delete [] m_data;
        m_data= new X[i*j];
    }
    /** return true if array is empty*/
    bool Empty() const
    {
        return ( m_data == NULL );
    }
    /** resize an existing array to i rows and j columns, erase all previous data and set all values to x*/
    void Initialize ( size_t i, size_t j,const X& x );
    /** destrctor, erase all elements from the heap*/
    ~D2Array()
    {
        if ( m_data )
            delete [] m_data;
    }
    void Clear()
    {
        if ( m_data )
            delete [] m_data;
        m_data=NULL;
        m_rows=m_columns=0;
    }
    /** @return reference to element at (i,j) position*/
    X& operator() ( size_t i, size_t j );
    /** @return const reference to element at (i,j) position*/
    const X& operator() ( size_t i, size_t j ) const;
    /** casting operator to X*/
    operator X*()
    {
        return m_data;
    }
    /** casting operator to X*/
    operator const X*() const
    {
        return m_data;
    }
    //arithmetic operators
    /** @return the number of rows*/
    size_t NRows() const
    {
        return m_rows;
    }
    /** @return the number of columns*/
    size_t NColumns() const
    {
        return m_columns;
    }
    /** @return the maximum value*/
    X Max() const;
    /** @return the minimum value*/
    X Min() const;
    /** swap two rows*/
    void SwapRows ( size_t i, size_t j );
    /** swap tow columns*/
    void SwapColumns ( size_t i, size_t j );
    /** remove a row*/
    D2Array<X> RemoveRow ( size_t row );
    /** remove a column*/
    D2Array<X> RemoveColumn ( size_t column );
    /** Get a row as a  new monodimensional array*/
    D1Array<X> GetRow ( size_t row ) const;
    /** Get a column as a new monodimensional array*/
    D1Array<X> GetColumn ( size_t column ) const;
    /** return the symmetric component of a square matrix*/
    D2Array<X> SymmetricPart() const;
    /** return the antisymmetric component of a square matrix*/
    D2Array<X> AntiSymmetricPart() const;
    /** transpose the matrix*/
    void Transpose();
    /** Get the trace of a square matrix*/
    X Trace() const;

private:
    size_t m_rows;
    size_t m_columns;
};



template<class X>
D2Array<X>::D2Array ( size_t i, size_t j, const X& x ) : m_rows ( i ), m_columns ( j )
{
    size_t total=i*j;
    m_data = new X[total];
    for ( size_t k=0;k<total;++k )
        m_data[k]=x;
}

template<class X>
void D2Array<X>::Initialize ( size_t i, size_t j, const X& x )
{
    if ( m_data )
        delete [] m_data;
    m_rows=i;
    m_columns=j;
    size_t total=i*j;
    m_data = new X[total];
    for ( size_t k=0;k<total;++k )
        m_data[k]=x;

}

template <class X>
D2Array<X>::D2Array ( const D2Array<X>& x )
{
    m_rows=x.NRows();
    m_columns=x.NColumns();
    if ( x.m_data != NULL )
    {
        m_data= new X[m_rows*m_columns];
        memcpy ( m_data,x,m_rows*m_columns*sizeof ( X ) );
    }
    else m_data = NULL;
}

template <class X>
D2Array<X>& D2Array<X>::operator= ( const D2Array& x )
{
    if ( m_data )
        delete [] m_data;
    m_rows=x.NRows();
    m_columns=x.NColumns();
    if ( x.m_data == NULL )
    {
        m_data=NULL;
    }
    else
    {
        m_data= new X[m_rows*m_columns];
        memcpy ( m_data,x,m_rows*m_columns*sizeof ( X ) );
    }
    return *this;
}

template <class X> void D2Array<X>::operator += ( const D2Array& x)
{
    assert ( this->m_rows == x.m_rows && this->m_columns == x.m_columns );
    for(size_t i=0;i<m_rows;++i)
        for(size_t j=0;j<m_columns;++j)
            (*this)(i,j)+=x(i,j);
}

template <class X> void D2Array<X>::operator -= ( const D2Array& x)
{
    assert ( this->m_rows == x.m_rows && this->m_columns == x.m_columns );
    for(size_t i=0;i<m_rows;++i)
        for(size_t j=0;j<m_columns;++j)
            (*this)(i,j)-=x(i,j);
}

template <class X> void D2Array<X>::operator *=(const X& x)
{
    for(size_t i=0;i<this->NRows();++i)
        for(size_t j=0;j<this->NColumns();++j)
        {
            (*this)(i,j)*=x;
        }
}

template <class X>
X& D2Array<X>::operator() ( size_t i, size_t j )
{
#ifdef ROW_MAJOR
    return m_data[i*m_columns+j];
#else

    return m_data[j*m_rows+i];
#endif

}

template <class X>
const X& D2Array<X>::operator() ( size_t i, size_t j ) const
{

#ifdef ROW_MAJOR
    return m_data[i*m_columns+j];
#else

    return m_data[j*m_rows+i];
#endif


}

template <class X>
X D2Array<X>::Max() const
{
    X max=(*this)(0,0);
    for (size_t i=0; i<m_rows; ++i)
        for (size_t j=0; j<m_columns; ++j)
            if ((*this)(i,j)>max)
                max = (*this)(i,j);

    return max;
}


template <class X>
X D2Array<X>::Min() const
{
    X min=(*this)(0,0);
    for (size_t i=0; i<m_rows; ++i)
        for (size_t j=0; j<m_columns; ++j)
            if ((*this)(i,j)<min)
                min = (*this)(i,j);

    return min;
}


template <class X>
void D2Array<X>::SwapRows ( size_t i, size_t j )
{
#ifdef ROW_MAJOR
    for ( size_t k=0;k<m_columns;k++ )
    {
        X tmp=m_data[m_columns*j+k];
        m_data[m_columns*j+k]=m_data[m_columns*i+k];
        m_data[m_columns*i+k]=tmp;
    }
#else
    for ( size_t k=0;k<m_columns;++k )
    {
        X tmp=m_data[k*m_rows+j];
        m_data[k*m_rows+j]=m_data[k*m_rows+i];
        m_data[k*m_rows+i]=tmp;
    }
#endif
}

template <class X>
void D2Array<X>::SwapColumns ( size_t i, size_t j )
{
    for ( size_t k=0;k<m_rows;++k )
    {
        X tmp= ( *this ) ( k,i );
        ( *this ) ( k,i ) = ( *this ) ( k,j );
        ( *this ) ( k,j ) =tmp;
    }
}


template <class X>
D2Array<X> D2Array<X>::RemoveRow ( size_t row )
{
    D2Array<X> nmat ( m_rows-1,m_columns );
    size_t k=0;
    for ( size_t i=0;i<m_rows;i++ )
    {
        if ( i != row )
        {
            for ( size_t j=0;j<m_columns;j++ )
                nmat ( k,j ) = ( *this ) ( i,j );
            k++;
        }

    }

    return nmat;
}


template <class X>
D2Array<X> D2Array<X>::RemoveColumn ( size_t column )
{
    D2Array<X> nmat ( m_rows,m_columns-1 );

    for ( size_t i=0;i<m_rows;i++ )
    {
        size_t k=0;
        for ( size_t j=0;j<m_columns;j++ )
        {
            if ( j!=column )
                nmat ( i,k++ ) = ( *this ) ( i,j );
        }
    }
    return nmat;
}

template <class X>
D1Array<X> D2Array<X>::GetRow ( size_t row ) const
{
    D1Array<X> fila ( m_columns );
    for ( size_t i=0;i<m_columns;i++ )
        fila ( i ) = ( *this ) ( row,i );

    return fila;
}

template <class X>
D1Array<X> D2Array<X>::GetColumn ( size_t column ) const
{
    D1Array<X> columna ( m_rows );
    for ( size_t i=0;i<m_rows;i++ )
        columna ( i ) = ( *this ) ( i,column );

    return columna;
}

template <class X>
D2Array<X> D2Array<X>::SymmetricPart() const
{
    D2Array<X> t(*this);
    for(size_t i=0;i<this->NRows();++i)
        for(size_t j=i+1;j<this->NColumns();++j)
        {
            X mean=(t(i,j)+t(j,i))/2;
            t(i,j)=t(j,i)=mean;
        }
    return t;
}

template <class X>
D2Array<X> D2Array<X>::AntiSymmetricPart() const
{
    D2Array<X> t(*this);
    for(size_t i=0;i<this->NRows();++i)
        for(size_t j=i;j<this->NColumns();++j)
        {
            if ( i == j )
                t(i,i)=0;
            else
            {
                X mean=(t(i,j)-t(j,i))/2;
                t(i,j)=mean;
                t(j,i)=-mean;
            }
        }
    return t;
}


template <class X>
void D2Array<X>::Transpose()
{
    D2Array<X>& a=(*this);
    assert ( a.NRows() == a.NColumns() );
    for(size_t i=0;i<a.NRows();++i)
        for(size_t j=0;j<a.NColumns();++j)
        {
            X tmp=a(i,j);
            a(i,j)=a(j,i);
            a(j,i)=tmp;
        }
}

template <class X>
X D2Array<X>::Trace() const
{
    assert ( this->NRows() == this->NColumns() );
    X trace=0;
    for(size_t i=0;i<this->NRows();++i)
    {
        trace+=(*this)(i,i);
    }
    return trace;
}


/** specializations*/
inline D2Array<double> operator * ( const D2Array<double>& a,const D2Array<double>& b )
{
    assert ( a.NColumns() == b.NRows() );
    D2Array<double> c ( a.NRows(),b.NColumns() );

#ifdef WITH_MKL 
    int stride;
    if ( a.NRows() > a.NColumns() ) stride= a.NRows();
    else stride=a.NColumns();
    cblas_dgemm ( CblasColMajor,CblasNoTrans,CblasNoTrans,a.NRows(),b.NColumns(),a.NRows(),1.0,a,stride,b,b.NColumns(),0.0,c,stride );
#endif

#ifdef WITH_VECLIB
    long int stride;
    if ( a.NRows() > a.NColumns() ) stride= a.NRows();
    else stride=a.NColumns();
    cblas_dgemm ( CblasColMajor,CblasNoTrans,CblasNoTrans,a.NRows(),b.NColumns(),a.NRows(),1.0,a,stride,b,b.NColumns(),0.0,c,stride );
#endif
#ifdef WITH_LAPACK

    char transa='N';
    char transb='N';
    int m=a.NRows();
    int n=b.NColumns();
    int k=a.NColumns();
    double alpha=1.0;
    int lda=a.NRows();
    int ldb=k;
    int ldc=m;
    double beta=0;

    dgemm(&transa,&transb,&m,&n,&k,&alpha,const_cast<double*> ( static_cast<const double*> ( a )),&lda,const_cast<double*> ( static_cast<const double*> ( b )),&ldb,&beta,c,&ldc);
#endif

    return c;

}

inline D1Array<double> operator * ( const D2Array<double>& a,const D1Array<double>& b )
{
    assert ( a.NRows() == b.size() );

    D1Array<double> c ( a.NRows() );

#ifdef WITH_MKL 
    cblas_dgemv ( CblasColMajor,CblasNoTrans,a.NRows(),a.NColumns(),1.0,a,a.NRows(),b,1,0.0,c,1 );
#endif

#ifdef WITH_VECLIB
    cblas_dgemv ( CblasColMajor,CblasNoTrans,a.NRows(),a.NColumns(),1.0,a,a.NRows(),b,1,0.0,c,1 );
#endif
#ifdef WITH_LAPACK
    char trans='N';
    int nrows=a.NRows();
    int ncolumns=a.NColumns();
    double alpha=1.0;
    int lda=a.NRows();
    int incx=1;
    double beta=0.0;
    int incy=1;
    dgemv(&trans,&nrows,&ncolumns,&alpha,const_cast<double*> ( static_cast<const double*> ( a ) ),&lda,const_cast<double*> ( static_cast<const double*> ( b ) ),&incx,&beta,c,&incy);
#endif

    return c;
}
#if defined WITH_MKL || defined WITH_VECLIB || defined WITH_LAPACK
/** diagonalize a symmetric matrix a, eigenvectors will be stored as a columns and eigenvalues will be stored in vector e*/
inline int Eigen ( D2Array<double>& a, D1Array<double>& e )
{
    char jobz='V';
    char uplo='U';
#ifdef WITH_VECLIB
    __CLPK_integer ldwork=a.NColumns() *a.NColumns();
    D1Array<double> work ( ldwork );
    __CLPK_integer info;
    __CLPK_integer n=a.NColumns();
    __CLPK_integer lda=a.NColumns();
    dsyev_ ( &jobz,&uplo,&n,a,&lda,e,work,&ldwork,&info );
#else
    int ldwork=a.NColumns() *a.NColumns();
    D1Array<double> work ( ldwork );
    int info;
    int n=a.NColumns();
    int lda=a.NColumns();
    dsyev ( &jobz,&uplo,&n,a,&lda,e,work,&ldwork,&info );
#endif
    return info;
}

inline void SVDFit(const D1Array<double>& res, const D2Array<double>& modelmatrix, D1Array<double>& solution)
{

#if defined WITH_MKL || defined WITH_LAPACK
    D2Array<double> A=modelmatrix;
    D1Array<double> b=res;
    int m=A.NRows();
    int n=A.NColumns();
    int nrhs=1;
    int lda=A.NRows();
    int ldb=res.size();
    D1Array<double> work ( 7*m );
    int ldwork=7*m;
    int info;
    D1Array<double> s ( A.NColumns() );
    int rank;
    double rcond=-1; //use machine precision
    dgelss ( &m,&n,&nrhs,A,&lda,b,&ldb,s,&rcond,&rank,work,&ldwork,&info );

    if ( info < 0 )
        throw kryomol::Exception("Illegal entry to dgelss routine in SVDFit");
    if ( info > 0 )
        throw kryomol::Exception("dgelss routine has not converged in SVDFit");

    for(size_t i=0;i<A.NColumns();++i)
        solution(i)=b(i);

    return;

#endif

#ifdef WITH_VECLIB
    D2Array<double> A=modelmatrix;
    D1Array<double> b=res;
    __CLPK_integer m=A.NRows();
    __CLPK_integer n=A.NColumns();
    __CLPK_integer nrhs=1;
    __CLPK_integer lda=A.NRows();
    __CLPK_integer ldb=b.size();
    D1Array<double> work ( 7*m );
    __CLPK_integer ldwork=7*m;
    __CLPK_integer info;
    D1Array<double> s ( A.NColumns() );
    __CLPK_integer rank;
    double rcond=-1; //use machine precision
    dgelss_ ( &m,&n,&nrhs,A,&lda,b,&ldb,s,&rcond,&rank,work,&ldwork,&info );
    if ( info < 0 )
        throw kryomol::Exception("Illegal entry to dgelss routine in SVDFit");
    if ( info > 0 )
        throw kryomol::Exception("dgelss routine has not converged in SVDFit");

    for(size_t i=0;i<A.NColumns();++i)
        solution(i)=b(i);

    return;
#endif
    throw kryomol::Exception("SVD fit not implemented");
}

inline int Invert ( D2Array<double>& a )
{

#ifdef WITH_VECLIB
    __CLPK_integer m=a.NRows();
    __CLPK_integer n=a.NColumns();
    __CLPK_integer lda=m;
    __CLPK_integer info;
    __CLPK_integer* ipiv = new __CLPK_integer[m];
    dgetrf_ ( &m,&n,a,&lda,ipiv,&info );
#else
    int m=a.NRows();
    int n=a.NColumns();
    int lda=m;
    int info;
    int* ipiv = new int[m];
    dgetrf ( &m,&n,a,&lda,ipiv,&info );
#endif
    D1Array<double> work ( n );
#ifdef WITH_VECLIB
    __CLPK_integer lwork=n;
    dgetri_ ( &n,a,&lda,ipiv,work,&lwork,&info );
#else
    int lwork=n;
    dgetri ( &n,a,&lda,ipiv,work,&lwork,&info );
#endif
    delete [] ipiv;
    std::cout << "info from dgetri " << info << std::endl;
    return info;

}


#endif

/** @param library name of the mathematical library (ccmath or Intel mkl) */
/** @param storage return the matrix storage system as "row_major" or "column_major"*/
TOOLS_API void  MathInfo ( std::string& library, std::string& storage );

template <class X> class  Scalar
{
public:
    Scalar() :  m_value(NULL) {}
    /** */
    ~Scalar() { delete m_value; }
    X Value() const { return *m_value; }

    /** Build error object with value v*/
    explicit Scalar(X v)
    {
        m_value =new X(v);
    }
    /** copy constructor*/
    Scalar(const Scalar& v)
    {

        if (  v )
        {
            m_value = new X(v.Value());

        }
        else m_value=NULL;
        //   if (!m_value && !e)  //Do nothing


    }
    /** assignment operator*/
    Scalar& operator = (const Scalar& e)
    {
        if ( &e != this)
        {
            if ( m_value && e )
            {
                *m_value=e.Value();

            }
            if ( m_value && !e)
            {
                delete m_value;
                m_value=NULL;
            }
            if ( !m_value && e)
            {
                m_value = new X(e.Value());

            }
            //   if (!m_value && !e)  //Do nothing

        }
        return *this;
    }


    /** asignment operator*/
    Scalar& operator = (X v)
    {
        {
            if ( m_value ) *m_value=v;
            else
                m_value = new X(v);
        }
        return *this;
    }
private:
    typedef void (Scalar::*bool_type) () const;
    void helperfunction() const {}
private:
    //static unity m_unity;
    X* m_value;
public:
    /** @return true if error has been defined*/
    operator bool_type () const
    {
        if ( m_value != 0 ) return &Scalar::helperfunction;
        else return 0;
    }
};



/** @brief three-dimensional array

     abstraction of memory-contiguos three-dimensional arrays
     */


template<class X>
class D3Array
{
public:
    /** Build an empty N-dimensional array */
    D3Array()
    {
        m_rows=m_columns=m_files=0;
        m_data=NULL;
    }
    /** Build an array with  i rows and j columns*/
    explicit D3Array ( size_t i, size_t j, size_t k ) : m_rows ( i ), m_columns ( j ), m_files ( k )
    {
        assert( i > 0 && j > 0 && k > 0 );
#ifdef WITH_SIMD
        if (( typeid(*m_data)==typeid(float) ) && ( (i%4==0) && (j%4==0) && (k%4==0) ))
        {
            m_data = (float*) _mm_malloc(i*j*k*sizeof(float),16);
        }
        else
            m_data = new X[i*j*k];
#else
        m_data = new X[i*j*k];
#endif
    }
    /** Build an array with  i rows and j columns*, with values x*/
    D3Array ( size_t i,size_t j, size_t k, const X& x );
    /** copy constructor*/
    D3Array ( const D3Array& x );
    /** assignment operator*/
    D3Array& operator= ( const D3Array& x );
    /** sum operator */
    void operator+=(const D3Array& x);
    /** difference operator*/
    void operator-=(const D3Array& x);
    /** multiply by scalar x*/
    void operator*=(const X& x);
    /** resize an existing array to i rows and j columns, erase all previous data*/
    void Initialize ( size_t i, size_t j, size_t k)
    {
        assert( i > 0 && j > 0 && k > 0 );
#ifdef WITH_SIMD
        if (( typeid(*m_data)==typeid(float) ) && ( (i%4==0) && (j%4==0) && (k%4==0) ))
        {
            m_rows=i;
            m_columns=j;
            m_files=k;
            if ( m_data )
                _mm_free(m_data);
            m_data= (float*) (_mm_malloc(i*j*k*sizeof(float),16));

        }
        else
        {
            m_rows=i;
            m_columns=j;
            m_files=k;
            if ( m_data )
                delete [] m_data;
            m_data= new X[i*j*k];
        }
#else
        m_rows=i;
        m_columns=j;
        m_files=k;
        if ( m_data )
            delete [] m_data;
        m_data= new X[i*j*k];
#endif

    }
    /** return true if array is empty*/
    bool Empty() const
    {
        return ( m_data == NULL );
    }
    /** resize an existing array to i rows and j columns, erase all previous data and set all values to x*/
    void Initialize ( size_t i, size_t j, size_t k, const X& x );
    /** destructor, erase all elements from the heap*/
    ~D3Array()
    {
#ifdef WITH_SIMD
        if (( typeid(*m_data)==typeid(float) ) && ( (m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) ))
        {
            if ( m_data )
                _mm_free(m_data);
        }
        else
        {
            if ( m_data )
                delete [] m_data;
        }
#else
        if ( m_data )
            delete [] m_data;
#endif

    }
    void Clear()
    {
#ifdef WITH_SIMD
        if (( typeid(*m_data)==typeid(float) ) && ( (m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) ))
        {
            if ( m_data )
                _mm_free(m_data);
            m_data=NULL;
            m_rows=m_columns=m_files=0;
        }
        else
        {
            if ( m_data )
                delete [] m_data;
            m_data=NULL;
            m_rows=m_columns=m_files=0;
        }
#else
        if ( m_data )
            delete [] m_data;
        m_data=NULL;
        m_rows=m_columns=m_files=0;
#endif
    }
    /** @return reference to element at (i,j) position*/
    X& operator() ( size_t i, size_t j, size_t k );
    /** @return const referpartial ence to element at (i,j) position*/
    const X& operator() ( size_t i, size_t j, size_t k ) const;
    /** casting operator to X*/
    operator X*()
    {
        return m_data;
    }
    /** casting operator to X*/
    operator const X*() const
    {
        return m_data;
    }
    //arithmetic operators
    /** @return the number of elements in dimension X*/
    size_t NX() const
    {
        return m_rows;
    }
    /** @return the number of elements in dimension Y*/
    size_t NY() const
    {
        return m_columns;
    }
    /** @return the number of elements in dimension Z*/
    size_t NZ() const
    {
        return m_files;
    }
    /** @return the maximum value of the matrix*/
    X Max() const;
    /** @return the minimum value of the matrix*/
    X Min() const;
    /** @return the matrix 2D for a value k in dimension Z*/
    D2Array<X> XY(size_t k) const;
    /** @return the matrix 2D for a value j in dimension Y*/
    D2Array<X> XZ(size_t j) const;
    /** @return the matrix 2D for a value i in dimension X*/
    D2Array<X> YZ(size_t i) const;
    /** set to zero all elements of the array*/
    void SetToZero ();
    /** calculate the Hamard product of two arrays*/
    void Hamard ( const D3Array& x, const D3Array& y );
    /** calculate the square difference of two arrays and sum to the current array*/
    void SquareDifference ( const X a, const D3Array& x, const D3Array& y );
    //                        /** calculate the array of exponentails used for calculating the atomic orbitals*/
    //                        void CalculateExponential(const X Nx, const X Ny, const X Nz, const X step, const X cx, const X cy, const X cz, const X alpha, const X xs);
    //                        /** calculate the orbital type S*/
    //                        void CalculateOrbitalS (  const D3Array& exponential);
    //                        /** calculate the orbital type P*/
    //                        void CalculateOrbitalPx (  const D3Array& exponential, const X cx, const X Nx, const X step);
    //                        void CalculateOrbitalPy (  const D3Array& exponential, const X cy, const X Ny, const X step);
    //                        void CalculateOrbitalPz (  const D3Array& exponential, const X cz, const X Nz, const X step);
    //                        /** calculate the orbital type D*/
    //                        void CalculateOrbitalDxx (  const D3Array& exponential, const X cx, const X Nx, const X step);
    //                        void CalculateOrbitalDyy (  const D3Array& exponential, const X cy, const X Ny, const X step);
    //                        void CalculateOrbitalDzz (  const D3Array& exponential, const X cz, const X Nz, const X step);
    //                        void CalculateOrbitalDxy (  const D3Array& exponential, const X cx, const X cy, const X Nx, const X Ny, const X step);
    //                        void CalculateOrbitalDxz (  const D3Array& exponential, const X cx, const X cz, const X Nx, const X Nz, const X step);
    //                        void CalculateOrbitalDyz (  const D3Array& exponential, const X cy, const X cz, const X Ny, const X Nz, const X step);
    //                        /** calculate the orbital type F*/
    //                        void CalculateOrbitalFxxx (  const D3Array& exponential, const X cx, const X Nx, const X step);
    //                        void CalculateOrbitalFyyy (  const D3Array& exponential, const X cy, const X Ny, const X step);
    //                        void CalculateOrbitalFzzz (  const D3Array& exponential, const X cz, const X Nz, const X step);
    //                        void CalculateOrbitalFxxy (  const D3Array& exponential, const X cx, const X cy, const X Nx, const X Ny, const X step);
    //                        void CalculateOrbitalFxyy (  const D3Array& exponential, const X cx, const X cy, const X Nx, const X Ny, const X step);
    //                        void CalculateOrbitalFxxz (  const D3Array& exponential, const X cx, const X cz, const X Nx, const X Nz, const X step);
    //                        void CalculateOrbitalFxzz (  const D3Array& exponential, const X cx, const X cz, const X Nx, const X Nz, const X step);
    //                        void CalculateOrbitalFyyz (  const D3Array& exponential, const X cy, const X cz, const X Ny, const X Nz, const X step);
    //                        void CalculateOrbitalFyzz (  const D3Array& exponential, const X cy, const X cz, const X Ny, const X Nz, const X step);
    //                        void CalculateOrbitalFxyz (  const D3Array& exponential, const X cx, const X cy, const X cz, const X Nx, const X Ny, const X Nz, const X step);
    //                        /** calculate the orbital type D with spherical harmonics*/
    //                        void CalculateOrbitalDY0 (  const D3Array& exponential, const X cx, const X cy, const X cz, const X Nx, const X Ny, const X Nz, const X step);
    //                        void CalculateOrbitalDY1 (  const D3Array& exponential, const X cx, const X cz, const X Nx, const X Nz, const X step);
    //                        void CalculateOrbitalDY2 (  const D3Array& exponential, const X cy, const X cz, const X Ny, const X Nz, const X step);
    //                        void CalculateOrbitalDY3 (  const D3Array& exponential, const X cx, const X cy, const X Nx, const X Ny, const X step);
    //                        void CalculateOrbitalDY4 (  const D3Array& exponential, const X cx, const X cy, const X Nx, const X Ny, const X step);
    //                        /** calculate the orbital type F with spherical harmonics*/
    //                        void CalculateOrbitalFY0 (  const D3Array& exponential, const X cx, const X cy, const X cz, const X Nx, const X Ny, const X Nz, const X step);
    //                        void CalculateOrbitalFY1 (  const D3Array& exponential, const X cx, const X cy, const X cz, const X Nx, const X Ny, const X Nz, const X step);
    //                        void CalculateOrbitalFY2 (  const D3Array& exponential, const X cx, const X cy, const X cz, const X Nx, const X Ny, const X Nz, const X step);
    //                        void CalculateOrbitalFY3 (  const D3Array& exponential, const X cx, const X cy, const X cz, const X Nx, const X Ny, const X Nz, const X step);
    //                        void CalculateOrbitalFY4 (  const D3Array& exponential, const X cx, const X cy, const X cz, const X Nx, const X Ny, const X Nz, const X step);
    //                        void CalculateOrbitalFY5 (  const D3Array& exponential, const X cx, const X cy, const X Nx, const X Ny, const X step);
    //                        void CalculateOrbitalFY6 (  const D3Array& exponential, const X cx, const X cy, const X Nx, const X Ny, const X step);

    //                        /** calculate the orbital type S*/
    //                        void CalculateOrbitalS ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        /** calculate the orbital type P*/
    //                        void CalculateOrbitalPx ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xp, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalPy ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xp, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalPz ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xp, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        /** calculate the orbital type D*/
    //                        void CalculateOrbitalDxx ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalDxy ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalDxz ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalDyy ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalDyz ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalDzz ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        /** calculate the orbital type F*/
    //                        void CalculateOrbitalFxxx ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalFxxy ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalFxxz ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalFxyy ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalFxyz ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalFxzz ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalFyyy ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalFyyz ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalFyzz ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalFzzz ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        /** calculate the orbital type D with spherical harmonics*/
    //                        void CalculateOrbitalDY0 ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalDY1 ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalDY2 ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalDY3 ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalDY4 ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        /** calculate the orbital type F with spherical harmonics*/
    //                        void CalculateOrbitalFY0 ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalFY1 ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalFY2 ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalFY3 ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalFY4 ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalFY5 ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);
    //                        void CalculateOrbitalFY6 ( const D3Array& x, const D3Array& y, const D3Array& z, const X& xs, const X& alpha, const X& cx, const X& cy, const X& cz);


protected:
    X* m_data;
    size_t m_rows;
    size_t m_columns;
    size_t m_files;
};



template<class X>
D3Array<X>::D3Array ( size_t i, size_t j, size_t h, const X& x ) : m_rows ( i ), m_columns ( j ), m_files ( h )
{
#ifdef WITH_SIMD
    if (( typeid(*m_data)==typeid(float) ) && ( (m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) ))
    {
        size_t total=i*j*h;
        m_data = (float*) (_mm_malloc(total*sizeof(float),16));
        __m128* ptr = (__m128*) (m_data);
        for (size_t k=0;k<total/4;++k)
            ptr[k] = _mm_set_ps1(x);
        m_data = (float*) ptr;
    }
    else
    {
        size_t total=i*j*h;
        m_data = new X[total];
        for ( size_t k=0;k<total;++k )
            m_data[k]=x;
    }
#else
    size_t total=i*j*h;
    m_data = new X[total];
    for ( size_t k=0;k<total;++k )
        m_data[k]=x;
#endif
}

template<class X>
void D3Array<X>::Initialize ( size_t i, size_t j, size_t h, const X& x )
{
#ifdef WITH_SIMD
    if (( typeid(*m_data)==typeid(float) ) && ( (i%4==0) && (j%4==0) && (h%4==0) ))
    {
        if ( m_data )
            _mm_free(m_data);
        m_rows=i;
        m_columns=j;
        m_files=h;
        size_t total=i*j*h;
        m_data = (float*) (_mm_malloc(total*sizeof(float),16));
        __m128* ptr = (__m128*) (m_data);
        for (size_t k=0;k<total/4;++k)
            ptr[k] = _mm_set_ps1(x);
        m_data = (float*) ptr;
    }
    else
    {
        if ( m_data )
            delete [] m_data;
        m_rows=i;
        m_columns=j;
        m_files=h;
        size_t total=i*j*h;
        m_data = new X[total];
        for ( size_t k=0;k<total;++k )
            m_data[k]=x;
    }
#else
    if ( m_data )
        delete [] m_data;
    m_rows=i;
    m_columns=j;
    m_files=h;
    size_t total=i*j*h;
    m_data = new X[total];
    for ( size_t k=0;k<total;++k )
        m_data[k]=x;
#endif

}

template <class X>
D3Array<X>::D3Array ( const D3Array<X>& x )
{
#ifdef WITH_SIMD
    if (( typeid(*m_data)==typeid(float) ) && ( (x.m_rows%4==0) && (x.m_columns%4==0) && (x.m_files%4==0) ))
    {
        m_rows=x.NX();
        m_columns=x.NY();
        m_files=x.NZ();
        size_t total=m_rows*m_columns*m_files;
        if ( x.m_data != NULL )
        {
            m_data= (float*) (_mm_malloc(total*sizeof(float),16));
            __m128* ptr = (__m128*) m_data;
            size_t i=0;
            for (size_t k=0;k<total/4;++k)
            {
                ptr[k] = _mm_set_ps(x.m_data[i+3],x.m_data[i+2],x.m_data[i+1],x.m_data[i]);
                i+=4;
            }
            m_data = (float*) (ptr);
        }
        else m_data = NULL;
    }
    else
    {
        m_rows=x.NX();
        m_columns=x.NY();
        m_files=x.NZ();
        size_t total=m_rows*m_columns*m_files;
        if ( x.m_data != NULL )
        {
            m_data= new X[total];
            memcpy ( m_data,x,m_rows*m_columns*m_files*sizeof ( X ) );
        }
        else m_data = NULL;
    }
#else
    m_rows=x.NX();
    m_columns=x.NY();
    m_files=x.NZ();
    size_t total=m_rows*m_columns*m_files;
    if ( x.m_data != NULL )
    {
        m_data= new X[total];
        memcpy ( m_data,x,m_rows*m_columns*m_files*sizeof ( X ) );
    }
    else m_data = NULL;
#endif
}


template <class X>
D3Array<X>& D3Array<X>::operator= ( const D3Array& x )
{
#ifdef WITH_SIMD
    if (( typeid(*m_data)==typeid(float) ) && ( (x.m_rows%4==0) && (x.m_columns%4==0) && (x.m_files%4==0) ))
    {
        if ( m_data )
            _mm_free(m_data);
        m_rows=x.NX();
        m_columns=x.NY();
        m_files=x.NZ();
        size_t total = m_rows*m_columns*m_files;
        if ( x.m_data == NULL )
        {
            m_data=NULL;
        }
        else
        {
            m_data= (float*) (_mm_malloc(m_rows*m_columns*m_files*sizeof(float),16));
            __m128* ptr = (__m128*) m_data;
            size_t i=0;
            for (size_t k=0;k<total/4;++k)
            {
                ptr[k] = _mm_set_ps(x.m_data[i+3],x.m_data[i+2],x.m_data[i+1],x.m_data[i]);
                i+=4;
            }
            m_data = (float*) (ptr);
        }
        return *this;

    }
    else
    {
        if ( m_data )
            delete [] m_data;
        m_rows=x.NX();
        m_columns=x.NY();
        m_files=x.NZ();
        if ( x.m_data == NULL )
        {
            m_data=NULL;
        }
        else
        {
            m_data= new X[m_rows*m_columns*m_files];
            memcpy ( m_data,x,m_rows*m_columns*m_files*sizeof ( X ) );
        }
        return *this;
    }

#else
    if ( m_data )
        delete [] m_data;
    m_rows=x.NX();
    m_columns=x.NY();
    m_files=x.NZ();
    if ( x.m_data == NULL )
    {
        m_data=NULL;
    }
    else
    {
        m_data= new X[m_rows*m_columns*m_files];
        memcpy ( m_data,x,m_rows*m_columns*m_files*sizeof ( X ) );
    }
    return *this;
#endif
}



template <class X> void D3Array<X>::operator += ( const D3Array& x)
{
    assert ( this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files);
#ifdef WITH_SIMD
    if (( typeid(*m_data)==typeid(float) ) && ( (x.m_rows%4==0) && (x.m_columns%4==0) && (x.m_files%4==0) ))
    {
        __m128* ptr = (__m128*) (m_data);
        __m128* x_ptr = (__m128*) (x.m_data);

        for (size_t i=0; i<(m_rows*m_columns*m_files)/4; ++i)
        {
            ptr[i] = _mm_add_ps(ptr[i],x_ptr[i]);
        }
        m_data = (float*) (ptr);

    }
    else
    {
        for(size_t i=0;i<m_rows;++i)
            for(size_t j=0;j<m_columns;++j)
                for(size_t k=0;k<m_files;++k)
                    (*this)(i,j,k)+=x(i,j,k);
    }

#else
    for(size_t i=0;i<m_rows;++i)
        for(size_t j=0;j<m_columns;++j)
            for(size_t k=0;k<m_files;++k)
                (*this)(i,j,k)+=x(i,j,k);
#endif
}


template <class X> void D3Array<X>::operator -= ( const D3Array& x)
{
    assert ( this->m_rows == x.m_rows && this->m_columns == x.m_columns && this->m_files == x.m_files);
#ifdef WITH_SIMD
    if (( typeid(*m_data)==typeid(float) ) && ( (x.m_rows%4==0) && (x.m_columns%4==0) && (x.m_files%4==0) ))
    {
        __m128* ptr = (__m128*) (m_data);
        __m128* x_ptr = (__m128*) (x.m_data);

        for (size_t i=0; i<(m_rows*m_columns*m_files)/4; ++i)
        {
            ptr[i] = _mm_sub_ps(ptr[i],x_ptr[i]);
        }
        m_data = (float*) (ptr);

    }
    else
    {
        for(size_t i=0;i<m_rows;++i)
            for(size_t j=0;j<m_columns;++j)
                for(size_t k=0;k<m_files;++k)
                    (*this)(i,j,k)-=x(i,j,k);
    }

#else
    for(size_t i=0;i<m_rows;++i)
        for(size_t j=0;j<m_columns;++j)
            for(size_t k=0;k<m_files;++k)
                (*this)(i,j,k)-=x(i,j,k);
#endif
}

template <class X> void D3Array<X>::operator *=(const X& x)
{
#ifdef WITH_SIMD
    if (( typeid(*m_data)==typeid(float) ) && ( (m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) ))
    {
        __m128* ptr = (__m128*) (m_data);
        __m128* x4 = (__m128*) _mm_malloc(4*sizeof(float),16);

        *x4 = _mm_set_ps1(x);
        for (size_t i=0; i<(m_rows*m_columns*m_files)/4; ++i)
        {
            ptr[i] = _mm_mul_ps(*x4,ptr[i]);
        }
        m_data = (float*) (ptr);

    }
    else
    {
        for(size_t i=0;i<m_rows;++i)
            for(size_t j=0;j<m_columns;++j)
                for(size_t k=0;k<m_files;++k)
                    (*this)(i,j,k)*=x;
    }

#else
    for(size_t i=0;i<m_rows;++i)
        for(size_t j=0;j<m_columns;++j)
            for(size_t k=0;k<m_files;++k)
                (*this)(i,j,k)*=x;
#endif
}

template <class X>
X& D3Array<X>::operator() ( size_t i, size_t j, size_t k )
{
#ifdef ROW_MAJOR
    return m_data[i*m_columns*m_files+j*m_files+k];
#else
    return m_data[k*m_rows*m_columns+j*m_rows+i];
#endif

}

template <class X>
X D3Array<X>::Max() const
{
    X max=(*this)(0,0,0);
    for (size_t i=0; i<m_rows; ++i)
        for (size_t j=0; j<m_columns; ++j)
            for (size_t k=0; k<m_columns; ++k)
                if ((*this)(i,j,k)>max)
                    max = (*this)(i,j,k);

    return max;
}

template <class X>
X D3Array<X>::Min() const
{
    X min=(*this)(0,0,0);
    for (size_t i=0; i<m_rows; ++i)
        for (size_t j=0; j<m_columns; ++j)
            for (size_t k=0; k<m_columns; ++k)
                if ((*this)(i,j,k)<min)
                    min = (*this)(i,j,k);

    return min;
}

template<class X>
void D3Array<X>::SetToZero()
{
#ifdef WITH_SIMD
    if (( typeid(*m_data)==typeid(float) ) && ( (m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) ))
    {
        __m128* ptr = (__m128*) (m_data);
        for (size_t k=0; k<m_rows*m_columns*m_files/4;++k)
            ptr[k] = _mm_set_ps1(0.0);
        m_data = (float*) ptr;
    }
    else
    {
        for ( size_t k=0;k<m_rows*m_columns*m_files;++k )
            m_data[k]=0;
    }
#else
    for ( size_t k=0;k<(m_rows*m_columns*m_files);++k )
        m_data[k]=0;
#endif
}

template<class X>
void D3Array<X>::Hamard ( const D3Array& x, const D3Array& y)
{
    assert ( this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files);
#ifdef WITH_SIMD
    if (( typeid(*m_data)==typeid(float) ) && ( (m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) ))
    {
        __m128* ptr = (__m128*) (m_data);
        __m128* x4 = (__m128*) (x.m_data);
        __m128* y4 = (__m128*) (y.m_data);
        for (size_t k=0; k<m_rows*m_columns*m_files/4;++k)
        {
            x4[k] = _mm_mul_ps(x4[k],y4[k]);
            ptr[k] = _mm_add_ps(ptr[k],x4[k]);
        }
        m_data = (float*) ptr;
    }
    else
    {
        for(size_t i=0;i<m_rows;++i)
            for(size_t j=0;j<m_columns;++j)
                for(size_t k=0;k<m_files;++k)
                    (*this)(i,j,k)+=x(i,j,k)*y(i,j,k);
    }
#else
    for(size_t i=0;i<m_rows;++i)
        for(size_t j=0;j<m_columns;++j)
            for(size_t k=0;k<m_files;++k)
                (*this)(i,j,k)+=x(i,j,k)*y(i,j,k);
#endif
}

template<class X>
void D3Array<X>::SquareDifference ( const X a, const D3Array& x, const D3Array& y)
{
    assert ( this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files);

#ifdef WITH_SIMD
    if (( typeid(*m_data)==typeid(float) ) && ( (m_rows%4==0) && (m_columns%4==0) && (m_files%4==0) ))
    {
        __m128* ptr = (__m128*) (m_data);
        __m128* x4 = (__m128*) (x.m_data);
        __m128* y4 = (__m128*) (y.m_data);
        __m128* sub4 = (__m128*) _mm_malloc(4*sizeof(float),16);
        __m128* a4 = (__m128*) _mm_malloc(4*sizeof(float),16);
        *a4 = _mm_set_ps1(a);
        for (size_t k=0; k<m_rows*m_columns*m_files/4;++k)
        {
            x4[k] = _mm_mul_ps(x4[k],x4[k]);
            y4[k] = _mm_mul_ps(y4[k],y4[k]);
            *sub4 = _mm_sub_ps(x4[k],y4[k]);
            *sub4 = _mm_mul_ps(*sub4,*a4);
            ptr[k] = _mm_add_ps(ptr[k],*sub4);
        }
        m_data = (float*) ptr;
    }
    else
    {
        for(size_t i=0;i<m_rows;++i)
            for(size_t j=0;j<m_columns;++j)
                for(size_t k=0;k<m_files;++k)
                    (*this)(i,j,k) += a*(x(i,j,k)*x(i,j,k) - y(i,j,k)*y(i,j,k));
    }
#else
    for(size_t i=0;i<m_rows;++i)
        for(size_t j=0;j<m_columns;++j)
            for(size_t k=0;k<m_files;++k)
                (*this)(i,j,k) += a*(x(i,j,k)*x(i,j,k) - y(i,j,k)*y(i,j,k));
#endif
}

template<class X>
D2Array<X> D3Array<X>::XY(size_t k) const
{
    D2Array<X> xy(m_rows,m_columns);
    if (k < m_files)
    {
        for (size_t i=0; i<m_rows; ++i)
            for (size_t j=0; j<m_columns; ++j)
                xy(i,j) = (*this)(i,j,k);
    }
    return xy;
}

template<class X>
D2Array<X> D3Array<X>::XZ(size_t j) const
{
    D2Array<X> xz(m_rows,m_files);
    if (j < m_columns)
    {
        for (size_t i=0; i<m_rows; ++i)
            for (size_t k=0; k<m_files; ++k)
                xz(i,k) = (*this)(i,j,k);
    }
    return xz;
}

template<class X>
D2Array<X> D3Array<X>::YZ(size_t i) const
{
    D2Array<X> yz(m_columns,m_files);
    if (i < m_rows)
    {
        for (size_t j=0; j<m_columns; ++j)
            for (size_t k=0; k<m_files; ++k)
                yz(j,k) = (*this)(i,j,k);
    }
    return yz;
}


template <class X>
const X& D3Array<X>::operator() ( size_t i, size_t j, size_t k ) const
{

#ifdef ROW_MAJOR
    return m_data[i*m_columns*m_files+j*m_files+k];
#else

    return m_data[k*m_rows*m_columns+j*m_rows+i];
#endif

}



/*       template<class X>
        void D3Array<X>::Hamard ( const D3Array<X>& x, const D3Array<X>& y)
        {
            assert ( this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files);
#ifdef WITH_CUDA
                if ( typeid(m_data).name() == "Pf" )
                {
                    scudahamard(x.m_data,y.m_data,m_data,m_data,m_rows*m_columns*m_files);

                }
                else
                {
                    for(size_t i=0;i<m_rows;++i)
                            for(size_t j=0;j<m_columns;++j)
                                for(size_t k=0;k<m_files;++k)
                                    (*this)(i,j,k)+=x(i,j,k)*y(i,j,k);
                }

#else

            for(size_t i=0;i<m_rows;++i)
                    for(size_t j=0;j<m_columns;++j)
                        for(size_t k=0;k<m_files;++k)
                            (*this)(i,j,k)+=x(i,j,k)*y(i,j,k);
#endif
        }

        template<class X>
        void D3Array<X>::SquareDifference ( const X a, const D3Array<X>& x, const D3Array<X>& y)
        {
            assert ( this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files);
#ifdef WITH_CUDA
                if ( typeid(m_data).name() == "Pf" )
                {
                    scudasquarediff(a,x.m_data,y.m_data,m_data,m_rows*m_columns*m_files);

                }
                else
                {
                    for(size_t i=0;i<m_rows;++i)
                            for(size_t j=0;j<m_columns;++j)
                                for(size_t k=0;k<m_files;++k)
                                    (*this)(i,j,k) += a*(x(i,j,k)*x(i,j,k) - y(i,j,k)*y(i,j,k));
                }

#else
            for(size_t i=0;i<m_rows;++i)
                    for(size_t j=0;j<m_columns;++j)
                        for(size_t k=0;k<m_files;++k)
                            (*this)(i,j,k) += a*(x(i,j,k)*x(i,j,k) - y(i,j,k)*y(i,j,k));
#endif
        }

        template<class X>
        void D3Array<X>::CalculateExponential(const X Nx, const X Ny, const X Nz, const X step, const X cx, const X cy, const X cz, const X alpha, const X xs) //const X &cx, const X &cy, const X &cz)
        {
            qDebug() << "Exponential: " << exp(-alpha) << endl;

            X x = -Nx/2;
            for(size_t i=0;i<m_rows;++i)
            {
                X y = -Ny/2;
                for(size_t j=0;j<m_columns;++j)
                {
                    X z = -Nz/2;
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) += xs*(exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz))));
                        z += step;
                    }
                    y += step;
                }
                x += step;
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalS( const D3Array &exponential)
        {
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
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalPx( const D3Array &exponential, const X cx, const X Nx, const X step)
        {
            X x = -Nx/2;
            for(size_t i=0;i<m_rows;++i)
            {
                for(size_t j=0;j<m_columns;++j)
                {
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(x-cx)*PC::AngstromToBohr;
                    }
                }
                x+=step;
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalPy( const D3Array &exponential, const X cy, const X Ny, const X step)
        {
            for(size_t i=0;i<m_rows;++i)
            {
                X y = -Ny/2;
                for(size_t j=0;j<m_columns;++j)
                {
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(y-cy)*PC::AngstromToBohr;
                    }
                    y+=step;
                }
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalPz( const D3Array &exponential, const X cz, const X Nz, const X step)
        {
            for(size_t i=0;i<m_rows;++i)
            {
                for(size_t j=0;j<m_columns;++j)
                {
                    X z = -Nz/2;
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(z-cz)*PC::AngstromToBohr;
                        z+=step;
                    }
                }
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalDxx( const D3Array &exponential, const X cx, const X Nx, const X step)
        {
            X x = -Nx/2;
            for(size_t i=0;i<m_rows;++i)
            {
                for(size_t j=0;j<m_columns;++j)
                {
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(x-cx)*PC::AngstromToBohr*(x-cx)*PC::AngstromToBohr;
                    }
                }
                x+=step;
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalDyy( const D3Array &exponential, const X cy, const X Ny, const X step)
        {
            for(size_t i=0;i<m_rows;++i)
            {
                X y = -Ny/2;
                for(size_t j=0;j<m_columns;++j)
                {
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(y-cy)*PC::AngstromToBohr*(y-cy)*PC::AngstromToBohr;
                    }
                    y+=step;
                }
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalDzz( const D3Array &exponential, const X cz, const X Nz, const X step)
        {
            for(size_t i=0;i<m_rows;++i)
            {
                for(size_t j=0;j<m_columns;++j)
                {
                    X z = -Nz/2;
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(z-cz)*PC::AngstromToBohr*(z-cz)*PC::AngstromToBohr;
                        z+=step;
                    }
                }
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalDxy( const D3Array &exponential, const X cx, const X cy, const X Nx, const X Ny, const X step)
        {
            X x = -Nx/2;
            for(size_t i=0;i<m_rows;++i)
            {
                X y = -Ny/2;
                for(size_t j=0;j<m_columns;++j)
                {
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(x-cx)*PC::AngstromToBohr*(y-cy)*PC::AngstromToBohr;
                    }
                    y+=step;
                }
                x+=step;
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalDxz( const D3Array &exponential, const X cx, const X cz, const X Nx, const X Nz, const X step)
        {
            X x = -Nx/2;
            for(size_t i=0;i<m_rows;++i)
            {
                for(size_t j=0;j<m_columns;++j)
                {
                    X z = -Nz/2;
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(x-cx)*PC::AngstromToBohr*(z-cz)*PC::AngstromToBohr;
                        z+=step;
                    }
                }
                x+=step;
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalDyz( const D3Array &exponential, const X cy, const X cz, const X Ny, const X Nz, const X step)
        {
            for(size_t i=0;i<m_rows;++i)
            {
                X y = -Ny/2;
                for(size_t j=0;j<m_columns;++j)
                {
                    X z = -Nz/2;
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(y-cy)*PC::AngstromToBohr*(z-cz)*PC::AngstromToBohr;
                        z += step;
                    }
                    y += step;
                }
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalFxxx( const D3Array &exponential, const X cx, const X Nx, const X step)
        {
            X x = -Nx/2;
            for(size_t i=0;i<m_rows;++i)
            {
                for(size_t j=0;j<m_columns;++j)
                {
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(x-cx)*PC::AngstromToBohr*(x-cx)*PC::AngstromToBohr*(x-cx)*PC::AngstromToBohr;
                    }
                }
                x += step;
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalFyyy( const D3Array &exponential, const X cy, const X Ny, const X step)
        {
            for(size_t i=0;i<m_rows;++i)
            {
                X y = -Ny/2;
                for(size_t j=0;j<m_columns;++j)
                {
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(y-cy)*PC::AngstromToBohr*(y-cy)*PC::AngstromToBohr*(y-cy)*PC::AngstromToBohr;
                    }
                    y += step;
                }
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalFzzz( const D3Array &exponential, const X cz, const X Nz, const X step)
        {
            for(size_t i=0;i<m_rows;++i)
            {
                for(size_t j=0;j<m_columns;++j)
                {
                    X z = -Nz/2;
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(z-cz)*PC::AngstromToBohr*(z-cz)*PC::AngstromToBohr*(z-cz)*PC::AngstromToBohr;
                        z += step;
                    }
                }
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalFxxy( const D3Array &exponential, const X cx, const X cy, const X Nx, const X Ny, const X step)
        {
            X x = -Nx/2;
            for(size_t i=0;i<m_rows;++i)
            {
                X y = -Ny/2;
                for(size_t j=0;j<m_columns;++j)
                {
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(x-cx)*PC::AngstromToBohr*(y-cy)*PC::AngstromToBohr*(x-cx)*PC::AngstromToBohr;
                    }
                    y += step;
                }
                x += step;
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalFxyy( const D3Array &exponential, const X cx, const X cy, const X Nx, const X Ny, const X step)
        {
            X x = -Nx/2;
            for(size_t i=0;i<m_rows;++i)
            {
                X y = -Ny/2;
                for(size_t j=0;j<m_columns;++j)
                {
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(x-cx)*PC::AngstromToBohr*(y-cy)*PC::AngstromToBohr*(y-cy)*PC::AngstromToBohr;
                    }
                    y += step;
                }
                x += step;
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalFxxz( const D3Array &exponential, const X cx, const X cz, const X Nx, const X Nz, const X step)
        {
            X x = -Nx/2;
            for(size_t i=0;i<m_rows;++i)
            {
                for(size_t j=0;j<m_columns;++j)
                {
                    X z = -Nz/2;
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(x-cx)*PC::AngstromToBohr*(x-cx)*PC::AngstromToBohr*(z-cz)*PC::AngstromToBohr;
                        z += step;
                    }
                }
                x += step;
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalFxzz( const D3Array &exponential, const X cx, const X cz, const X Nx, const X Nz, const X step)
        {
            X x = -Nx/2;
            for(size_t i=0;i<m_rows;++i)
            {
                for(size_t j=0;j<m_columns;++j)
                {
                    X z = -Nz/2;
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(x-cx)*PC::AngstromToBohr*(z-cz)*PC::AngstromToBohr*(z-cz)*PC::AngstromToBohr;
                        z += step;
                    }
                }
                x += step;
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalFyyz( const D3Array &exponential, const X cy, const X cz, const X Ny, const X Nz, const X step)
        {
            for(size_t i=0;i<m_rows;++i)
            {
                X y = -Ny/2;
                for(size_t j=0;j<m_columns;++j)
                {
                    X z = -Nz/2;
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(y-cy)*PC::AngstromToBohr*(y-cy)*PC::AngstromToBohr*(z-cz)*PC::AngstromToBohr;
                        z += step;
                    }
                    y += step;
                }
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalFyzz( const D3Array &exponential, const X cy, const X cz, const X Ny, const X Nz, const X step)
        {
            for(size_t i=0;i<m_rows;++i)
            {
                X y = -Ny/2;
                for(size_t j=0;j<m_columns;++j)
                {
                    X z = -Nz/2;
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(z-cz)*PC::AngstromToBohr*(y-cy)*PC::AngstromToBohr*(z-cz)*PC::AngstromToBohr;
                        z += step;
                    }
                    y += step;
                }
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalFxyz( const D3Array &exponential, const X cx, const X cy, const X cz, const X Nx, const X Ny, const X Nz, const X step)
        {
            X x = -Nx/2;
            for(size_t i=0;i<m_rows;++i)
            {
                X y = -Ny/2;
                for(size_t j=0;j<m_columns;++j)
                {
                    X z = -Nz/2;
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(x-cx)*PC::AngstromToBohr*(y-cy)*PC::AngstromToBohr*(z-cz)*PC::AngstromToBohr;
                        z += step;
                    }
                    y += step;
                }
                x += step;
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalDY0( const D3Array &exponential, const X cx, const X cy, const X cz, const X Nx, const X Ny, const X Nz, const X step)
        {
            X m_pi = 3.14159265358979323846;
            X x = -Nx/2;
            for(size_t i=0;i<m_rows;++i)
            {
                X y = -Ny/2;
                for(size_t j=0;j<m_columns;++j)
                {
                    X z = -Nz/2;
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*0.25*sqrt(5/m_pi)*(-(x-cx)*(x-cx)-(y-cy)*(y-cy)+2*(z-cz)*(z-cz))*PC::AngstromToBohr*PC::AngstromToBohr;
                        z += step;
                    }
                    y += step;
                }
                x += step;
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalDY1( const D3Array &exponential, const X cx, const X cz, const X Nx, const X Nz, const X step)
        {
            X m_pi = 3.14159265358979323846;
            X x = -Nx/2;
            for(size_t i=0;i<m_rows;++i)
            {
                for(size_t j=0;j<m_columns;++j)
                {
                    X z = -Nz/2;
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(0.5)*sqrt(15/(m_pi))*(x-cx)*(z-cz)*PC::AngstromToBohr*PC::AngstromToBohr;
                        z += step;
                    }
                }
                x += step;
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalDY2( const D3Array &exponential, const X cy, const X cz, const X Ny, const X Nz, const X step)
        {
            X m_pi = 3.14159265358979323846;
            for(size_t i=0;i<m_rows;++i)
            {
                X y = -Ny/2;
                for(size_t j=0;j<m_columns;++j)
                {
                    X z = -Nz/2;
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(0.5)*sqrt(15/(m_pi))*(y-cy)*(z-cz)*PC::AngstromToBohr*PC::AngstromToBohr;
                        z += step;
                    }
                    y += step;
                }
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalDY3( const D3Array &exponential, const X cx, const X cy, const X Nx, const X Ny, const X step)
        {
            X m_pi = 3.14159265358979323846;
            X x = -Nx/2;
            for(size_t i=0;i<m_rows;++i)
            {
                X y = -Ny/2;
                for(size_t j=0;j<m_columns;++j)
                {
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*0.25*sqrt(15/(m_pi))*((x-cx)*(x-cx)-(y-cy)*(y-cy))*PC::AngstromToBohr*PC::AngstromToBohr;
                    }
                    y += step;
                }
                x += step;
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalDY4( const D3Array &exponential, const X cx, const X cy, const X Nx, const X Ny, const X step)
        {
            X m_pi = 3.14159265358979323846;
            X x = -Nx/2;
            for(size_t i=0;i<m_rows;++i)
            {
                X y = -Ny/2;
                for(size_t j=0;j<m_columns;++j)
                {
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*0.5*sqrt(15/(m_pi))*((x-cx)*(y-cy))*PC::AngstromToBohr*PC::AngstromToBohr;
                    }
                    y += step;
                }
                x += step;
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalFY0( const D3Array &exponential, const X cx, const X cy, const X cz, const X Nx, const X Ny, const X Nz, const X step)
        {
            X m_pi = 3.14159265358979323846;
            X x = -Nx/2;
            for(size_t i=0;i<m_rows;++i)
            {
                X y = -Ny/2;
                for(size_t j=0;j<m_columns;++j)
                {
                    X z = -Nz/2;
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*0.25*sqrt(7/m_pi)*(2*(z-cz)*(z-cz)-3*(x-cx)*(x-cx)-3*(y-cy)*(y-cy))*(z-cz)*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
                        z += step;
                    }
                    y += step;
                }
                x += step;
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalFY1( const D3Array &exponential, const X cx, const X cy, const X cz, const X Nx, const X Ny, const X Nz, const X step)
        {
            X m_pi = 3.14159265358979323846;
            X x = -Nx/2;
            for(size_t i=0;i<m_rows;++i)
            {
                X y = -Ny/2;
                for(size_t j=0;j<m_columns;++j)
                {
                    X z = -Nz/2;
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(0.25)*sqrt(21/2*m_pi)*(4*(z-cz)*(z-cz)-(x-cx)*(x-cx)-(y-cy)*(y-cy))*(x-cx)*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
                        z += step;
                    }
                    y += step;
                }
                x += step;
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalFY2( const D3Array &exponential, const X cx, const X cy, const X cz, const X Nx, const X Ny, const X Nz, const X step)
        {
            X m_pi = 3.14159265358979323846;
            X x = -Nx/2;
            for(size_t i=0;i<m_rows;++i)
            {
                X y = -Ny/2;
                for(size_t j=0;j<m_columns;++j)
                {
                    X z = -Nz/2;
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(0.25)*sqrt(21/2*m_pi)*(4*(z-cz)*(z-cz)-(x-cx)*(x-cx)-(y-cy)*(y-cy))*(y-cy)*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
                        z += step;
                    }
                    y += step;
                }
                x += step;
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalFY3( const D3Array &exponential, const X cx, const X cy, const X cz, const X Nx, const X Ny, const X Nz, const X step)
        {
            X m_pi = 3.14159265358979323846;
            X x = -Nx/2;
            for(size_t i=0;i<m_rows;++i)
            {
                X y = -Ny/2;
                for(size_t j=0;j<m_columns;++j)
                {
                    X z = -Nz/2;
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*0.25*sqrt(105/(m_pi))*((x-cx)*(x-cx)-(y-cy)*(y-cy))*(z-cz)*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
                        z += step;
                    }
                    y += step;
                }
                x += step;
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalFY4( const D3Array &exponential, const X cx, const X cy, const X cz, const X Nx, const X Ny, const X Nz, const X step)
        {
            X m_pi = 3.14159265358979323846;
            X x = -Nx/2;
            for(size_t i=0;i<m_rows;++i)
            {
                X y = -Ny/2;
                for(size_t j=0;j<m_columns;++j)
                {
                    X z = -Nz/2;
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*0.5*sqrt(105/(m_pi))*(x-cx)*(y-cy)*(z-cz)*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
                        z += step;
                    }
                    y += step;
                }
                x += step;
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalFY5( const D3Array &exponential, const X cx, const X cy, const X Nx, const X Ny, const X step)
        {
            X m_pi = 3.14159265358979323846;
            X x = -Nx/2;
            for(size_t i=0;i<m_rows;++i)
            {
                X y = -Ny/2;
                for(size_t j=0;j<m_columns;++j)
                {
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(0.25)*sqrt(35/2*m_pi)*((x-cx)*(x-cx)-3*(y-cy)*(y-cy))*(x-cx)*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
                    }
                    y += step;
                }
                x += step;
            }
        }

        template<class X>
        void D3Array<X>::CalculateOrbitalFY6( const D3Array &exponential, const X cx, const X cy, const X Nx, const X Ny, const X step)
        {
            X m_pi = 3.14159265358979323846;
            X x = -Nx/2;
            for(size_t i=0;i<m_rows;++i)
            {
                X y = -Ny/2;
                for(size_t j=0;j<m_columns;++j)
                {
                    for(size_t k=0;k<m_files;++k)
                    {
                        (*this)(i,j,k) = exponential(i,j,k)*(0.25)*sqrt(35/2*m_pi)*(3*(x-cx)*(x-cx)-(y-cy)*(y-cy))*(y-cy)*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
                    }
                    y += step;
                }
                x += step;
            }
        }*/


//------------------------------- Template class Node and Fifo for creating a FIFO stack ----------------------------//

template <class X, class Y>
class Node
{
public:
    Node() {}
    /** Build FIFO stack with a maximum of max elements*/
    explicit Node ( X x, Y y)
    {
        m_id = x;
        m_data = y;
    }
    /** copy constructor*/
    Node ( const Node& node );
    /** destructor, erase all elements from the heap*/
    ~Node()
    {}
    /** assignment operator*/
    Node& operator= ( const Node& x );
    /** @return ID*/
    X ID() const
    {
        return m_id;
    }
    /** @return Data*/
    Y Data() const
    {
        return m_data;
    }

private:
    X m_id;
    Y m_data;

};

template <class X, class Y>
Node<X,Y>::Node ( const Node& n )
{
    m_id = n.ID();
    m_data = n.Data();
}

template <class X, class Y>
Node<X,Y>& Node<X,Y>::operator= ( const Node& n )
{

    m_id = n.ID();
    m_data = n.Data();

    return *this;
}


/** @brief FIFO stack

     abstraction of memory-contiguos FIFO stack
     */

template <class X, class Y>
class Fifo
{
public:
    /** Build an empty FIFO stack */
    Fifo()
    {
        m_size = 0;
        m_nodes = NULL;
    }
    /** Build FIFO stack with a maximum of max elements*/
    explicit Fifo ( size_t max )
    {
        m_size = 0;
        m_nodes = NULL;
        m_max = max;
    }
    /** copy constructor*/
    //Fifo ( const Fifo& f );
    /** destructor, erase all elements from the heap*/
    ~Fifo()
    {
        if ( m_nodes != NULL )
            delete [] m_nodes;
    }
    /** assignment operator*/
    //Fifo& operator= ( const Fifo& f );
    /** return true if FIFO stack is empty*/
    bool Empty() const
    {
        return ( m_nodes == NULL );
    }

    void Clear()
    {
        if ( m_nodes )
            delete [] m_nodes;
        m_nodes=NULL;
        m_size=0;
    }
    /** @return reference to element with an ID id*/
    Node<X,Y>& operator() ( X id );
    /** @return the position in the list of the element with an ID id*/
    size_t find ( X id );
    /** push function */
    void push(const X x, const Y y);
    /** @return the number of elements in Fifo*/
    size_t Size() const
    {
        return m_size;
    }
    /** @return the maximum number of elements in FIFO stack*/
    size_t Max() const
    {
        return m_max;
    }
    /** set the maximum number of elements in FIFO stack*/
    void SetMax(size_t max)
    {
        m_max = max;
    }
    /** @return the list of Nodes of the FIFO stack*/
    Node<X,Y>* Nodes() const
    {
        return m_nodes;
    }


private:
    Node<X,Y>* m_nodes;
    size_t m_size;
    size_t m_max;
};

//template <class X, class Y>
//Fifo<X,Y>::Fifo ( const Fifo& f )
//{
//        m_size = f.Size();
//        m_max = f.Max();
//        if ( f.Nodes() != NULL )
//        {
//            m_nodes = new Node<X,Y>[m_size];
//            memcpy ( m_nodes,f,m_size*sizeof ( Node<X,Y> ) );
//        }
//        else m_nodes = NULL;
//}

//template <class X, class Y>
//Fifo<X,Y>& Fifo<X,Y>::operator= ( const Fifo& f )
//{
//        if ( m_nodes )
//                delete [] m_nodes;
//        m_size = f.Size();
//        if ( f.Nodes() == NULL )
//        {
//          m_nodes = NULL;
//        }
//        else
//        {
//            m_nodes= new Node<X,Y>[m_size];
//            memcpy ( m_nodes,f,m_size*sizeof ( Node<X,Y> ) );
//        }
//        return *this;
//}

template <class X, class Y>
Node<X,Y>& Fifo<X,Y>::operator() (X x)
{
    Node<X,Y> y;
    for (size_t it=0; it<m_size;++it)
    {
        Node<X,Y>& node = m_nodes[it];
        if (node.ID() == x)
            y = node;
    }
    return y;
}

template <class X, class Y>
size_t Fifo<X,Y>::find (X x)
{
    size_t pos = -1;
    for (size_t it=0; it<m_size; ++it)
    {
        Node<X,Y>& node = m_nodes[it];
        if (node.ID() == x)
            pos = it;
    }
    return pos;
}

template <class X, class Y>
void Fifo<X,Y>::push(const X x, const Y y)
{
    if (this->find(x) < 0)
    {
        if (m_size == m_max)
        {
            for (size_t i=m_size-1; i>0; --i)
            {
                m_nodes[i]=m_nodes[i-1];
            }
            m_nodes[0] = Node<X,Y>(x,y);
        }
        else
        {
            ++m_size;
            for (size_t i=m_size-1; i>0; --i)
            {
                m_nodes[i]=m_nodes[i-1];
            }
            m_nodes[0]=Node<X,Y>(x,y);
        }
    }
    else
    {
        size_t pos = this->find(x);
        if (pos!=0)
        {
            Node<X,Y> node = m_nodes[pos];
            for (size_t i=pos; i>0; --i)
            {
                m_nodes[i] = m_nodes[i-1];
            }
            m_nodes[0]=node;
        }
    }
}

}


//        template<class X>
//        void D3Array<X>::CalculateOrbitalS(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitals(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                    (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))));
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                            (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))));

//#endif
//        }

//        template<class X>
//        void D3Array<X>::CalculateOrbitalPx(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitalpx(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                    (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*(x(i,j,k)-cx);//*PC::AngstromToBohr;
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                            (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*(x(i,j,k)-cx);//*PC::AngstromToBohr;
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalPy(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitalpy(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                    (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*(y(i,j,k)-cy);//*PC::AngstromToBohr;
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                            (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*(y(i,j,k)-cy);//*PC::AngstromToBohr;
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalPz(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitalpz(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                    (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*(z(i,j,k)-cz);//*PC::AngstromToBohr;
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                            (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*(z(i,j,k)-cz);//*PC::AngstromToBohr;
//#endif
//        }

//        template<class X>
//        void D3Array<X>::CalculateOrbitalDxx(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitaldxx(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                    (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Tan30*(x(i,j,k)-cx)*(x(i,j,k)-cx);//*PC::AngstromToBohr*PC::AngstromToBohr;
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                            (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Tan30*(x(i,j,k)-cx)*(x(i,j,k)-cx);//*PC::AngstromToBohr*PC::AngstromToBohr;
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalDxy(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitaldxy(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                    (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*(x(i,j,k)-cx)*(y(i,j,k)-cy);//*PC::AngstromToBohr*PC::AngstromToBohr;
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                            (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*(x(i,j,k)-cx)*(y(i,j,k)-cy);//*PC::AngstromToBohr*PC::AngstromToBohr;
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalDxz(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitaldxz(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                    (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*(x(i,j,k)-cx)*(z(i,j,k)-cz);//*PC::AngstromToBohr*PC::AngstromToBohr;
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                            (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*(x(i,j,k)-cx)*(z(i,j,k)-cz);//*PC::AngstromToBohr*PC::AngstromToBohr;
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalDyy(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitaldyy(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                    (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Tan30*(y(i,j,k)-cy)*(y(i,j,k)-cy);//*PC::AngstromToBohr*PC::AngstromToBohr;
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                            (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Tan30*(y(i,j,k)-cy)*(y(i,j,k)-cy);//*PC::AngstromToBohr*PC::AngstromToBohr;
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalDyz(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitaldyz(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                    (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*(y(i,j,k)-cy)*(z(i,j,k)-cz);//*PC::AngstromToBohr*PC::AngstromToBohr;
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                            (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*(y(i,j,k)-cy)*(z(i,j,k)-cz);//*PC::AngstromToBohr*PC::AngstromToBohr;

//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalDzz(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitaldzz(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                    (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Tan30*(z(i,j,k)-cz)*(z(i,j,k)-cz);//*PC::AngstromToBohr*PC::AngstromToBohr;
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                            (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Tan30*(z(i,j,k)-cz)*(z(i,j,k)-cz);//*PC::AngstromToBohr*PC::AngstromToBohr;

//#endif
//        }

//        template<class X>
//        void D3Array<X>::CalculateOrbitalFxxx(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitalfxxx(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                    (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Sqrt15Div15*(x(i,j,k)-cx)*(x(i,j,k)-cx)*(x(i,j,k)-cx);//*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                            (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Sqrt15Div15*(x(i,j,k)-cx)*(x(i,j,k)-cx)*(x(i,j,k)-cx);//*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalFxxy(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitalfxxy(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                    (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Tan30*(x(i,j,k)-cx)*(x(i,j,k)-cx)*(y(i,j,k)-cy);//*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                            (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Tan30*(x(i,j,k)-cx)*(x(i,j,k)-cx)*(y(i,j,k)-cy);//*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalFxxz(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitalfxxz(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                    (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Tan30*(x(i,j,k)-cx)*(x(i,j,k)-cx)*(z(i,j,k)-cz);//*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                            (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Tan30*(x(i,j,k)-cx)*(x(i,j,k)-cx)*(z(i,j,k)-cz);//*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalFxyy(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitalfxyy(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                    (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Tan30*(x(i,j,k)-cx)*(y(i,j,k)-cy)*(y(i,j,k)-cy);//*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                            (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Tan30*(x(i,j,k)-cx)*(y(i,j,k)-cy)*(y(i,j,k)-cy);//*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalFxyz(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitalfxyz(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                    (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*(x(i,j,k)-cx)*(y(i,j,k)-cy)*(z(i,j,k)-cz);//*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                            (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*(x(i,j,k)-cx)*(y(i,j,k)-cy)*(z(i,j,k)-cz);//*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalFxzz(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitalfxzz(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                    (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Tan30*(x(i,j,k)-cx)*(z(i,j,k)-cz)*(z(i,j,k)-cz);//*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                            (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Tan30*(x(i,j,k)-cx)*(z(i,j,k)-cz)*(z(i,j,k)-cz);//*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalFyyy(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitalfyyy(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                    (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Sqrt15Div15*(y(i,j,k)-cy)*(y(i,j,k)-cy)*(y(i,j,k)-cy);//*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                            (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Sqrt15Div15*(y(i,j,k)-cy)*(y(i,j,k)-cy)*(y(i,j,k)-cy);//*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalFyyz(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitalfyyz(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                    (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Tan30*(y(i,j,k)-cy)*(y(i,j,k)-cy)*(z(i,j,k)-cz);//*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                            (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Tan30*(y(i,j,k)-cy)*(y(i,j,k)-cy)*(z(i,j,k)-cz);//*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalFyzz(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitalfyzz(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                    (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Tan30*(y(i,j,k)-cy)*(z(i,j,k)-cz)*(z(i,j,k)-cz);//*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                            (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Tan30*(y(i,j,k)-cy)*(z(i,j,k)-cz)*(z(i,j,k)-cz);//*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalFzzz(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitalfzzz(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                    (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Sqrt15Div15*(z(i,j,k)-cz)*(z(i,j,k)-cz)*(z(i,j,k)-cz);//*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                            (*this)(i,j,k) = (xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz))))*PC::Sqrt15Div15*(z(i,j,k)-cz)*(z(i,j,k)-cz)*(z(i,j,k)-cz);//*PC::AngstromToBohr*PC::AngstromToBohr*PC::AngstromToBohr;
//#endif
//        }

//        template<class X>
//        void D3Array<X>::CalculateOrbitalDY0(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//            X m_pi = 3.14159265358979323846;
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitaldy0(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                {
//                                    //X r2 = ((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz));
//                                    (*this)(i,j,k) = xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*0.25*sqrt(5/m_pi)*(-(x(i,j,k)-cx)*(x(i,j,k)-cx)-(y(i,j,k)-cy)*(y(i,j,k)-cy)+2*(z(i,j,k)-cz)*(z(i,j,k)-cz));
//                                }
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                        {
//                            //X r2 = (x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz);
//                            (*this)(i,j,k) = xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*0.25*sqrt(5/m_pi)*(-(x(i,j,k)-cx)*(x(i,j,k)-cx)-(y(i,j,k)-cy)*(y(i,j,k)-cy)+2*(z(i,j,k)-cz)*(z(i,j,k)-cz));
//                        }
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalDY1(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//            X m_pi = 3.14159265358979323846;
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitaldy1(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                {
//                                    //X r2 = (x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz);
//                                    (*this)(i,j,k) = xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*(-0.5)*sqrt(15/(2*m_pi))*(x(i,j,k)-cx)*(z(i,j,k)-cz);
//                                }
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                        {
//                            //X r2 = (x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz);
//                            (*this)(i,j,k) = (X) xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*(-0.5)*sqrt(15/(2*m_pi))*(x(i,j,k)-cx)*(z(i,j,k)-cz);
//                        }
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalDY2(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//            X m_pi = 3.14159265358979323846;
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitaldy2(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                {
//                                    //X r2 = (x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz);
//                                    (*this)(i,j,k) = (X) xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*(0.5)*sqrt(15/(2*m_pi))*(x(i,j,k)-cx)*(z(i,j,k)-cz);
//                                }
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                        {
//                            //X r2 = (x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz);
//                            (*this)(i,j,k) = (X) xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*(0.5)*sqrt(15/(2*m_pi))*(x(i,j,k)-cx)*(z(i,j,k)-cz);
//                        }
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalDY3(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//            X m_pi = 3.14159265358979323846;
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitaldy3(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                {
//                                    //X r2 = (x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz);
//                                    (*this)(i,j,k) = xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*0.25*sqrt(15/(2*m_pi))*((x(i,j,k)-cx)*(x(i,j,k)-cx)-(y(i,j,k)-cy)*(y(i,j,k)-cy));
//                                }
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                        {
//                            //X r2 = (x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz);
//                            (*this)(i,j,k) = xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*0.25*sqrt(15/(2*m_pi))*((x(i,j,k)-cx)*(x(i,j,k)-cx)-(y(i,j,k)-cy)*(y(i,j,k)-cy));
//                        }
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalDY4(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//            X m_pi = 3.14159265358979323846;
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitaldy4(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                {
//                                    //X r2 = (x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz);
//                                    (*this)(i,j,k) = xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*0.25*sqrt(15/(2*m_pi))*((x(i,j,k)-cx)*(x(i,j,k)-cx)-(y(i,j,k)-cy)*(y(i,j,k)-cy));
//                                }
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                        {
//                            //X r2 = (x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz);
//                            (*this)(i,j,k) = xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*0.25*sqrt(15/(2*m_pi))*((x(i,j,k)-cx)*(x(i,j,k)-cx)-(y(i,j,k)-cy)*(y(i,j,k)-cy));
//                        }
//#endif
//        }

//        template<class X>
//        void D3Array<X>::CalculateOrbitalFY0(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//            X m_pi = 3.14159265358979323846;
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitalfy0(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                {
//                                    //X r = sqrt((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz));
//                                    (*this)(i,j,k) = xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*0.25*sqrt(7/m_pi)*(2*(z(i,j,k)-cz)*(z(i,j,k)-cz)-3*(x(i,j,k)-cx)*(x(i,j,k)-cx)-3*(y(i,j,k)-cy)*(y(i,j,k)-cy))*(z(i,j,k)-cz)*(z(i,j,k)-cz);
//                                }
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                        {
//                            //X r = sqrt((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz));
//                            (*this)(i,j,k) = xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*0.25*sqrt(7/m_pi)*(2*(z(i,j,k)-cz)*(z(i,j,k)-cz)-3*(x(i,j,k)-cx)*(x(i,j,k)-cx)-3*(y(i,j,k)-cy)*(y(i,j,k)-cy))*(z(i,j,k)-cz)*(z(i,j,k)-cz);
//                        }
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalFY1(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//            X m_pi = 3.14159265358979323846;
//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitalfy1(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                {
//                                    //X r = sqrt((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz));
//                                    (*this)(i,j,k) = xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*(-0.125)*sqrt(21/m_pi)*(x(i,j,k)-cx)*(4*(z(i,j,k)-cz)*(z(i,j,k)-cz)-(x(i,j,k)-cx)*(x(i,j,k)-cx)-(y(i,j,k)-cy)*(y(i,j,k)-cy));
//                                }
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                        {
//                            //X r = sqrt((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz));
//                            (*this)(i,j,k) = xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*(-0.125)*sqrt(21/m_pi)*(x(i,j,k)-cx)*(4*(z(i,j,k)-cz)*(z(i,j,k)-cz)-(x(i,j,k)-cx)*(x(i,j,k)-cx)-(y(i,j,k)-cy)*(y(i,j,k)-cy));
//                        }
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalFY2(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//            X m_pi = 3.14159265358979323846;

//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitalfy2(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                {
//                                    //X r = sqrt((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz));
//                                    (*this)(i,j,k) = xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*(0.125)*sqrt(21/m_pi)*(x(i,j,k)-cx)*(4*(z(i,j,k)-cz)*(z(i,j,k)-cz)-(x(i,j,k)-cx)*(x(i,j,k)-cx)-(y(i,j,k)-cy)*(y(i,j,k)-cy));
//                                }
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                        {
//                            //X r = sqrt((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz));
//                            (*this)(i,j,k) = xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*(0.125)*sqrt(21/m_pi)*(x(i,j,k)-cx)*(4*(z(i,j,k)-cz)*(z(i,j,k)-cz)-(x(i,j,k)-cx)*(x(i,j,k)-cx)-(y(i,j,k)-cy)*(y(i,j,k)-cy));
//                        }
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalFY3(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//            X m_pi = 3.14159265358979323846;

//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitalfy3(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                {
//                                    //X r = sqrt((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz));
//                                    (*this)(i,j,k) = xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*0.25*sqrt(105/(2*m_pi))*((x(i,j,k)-cx)*(x(i,j,k)-cx)-(y(i,j,k)-cy)*(y(i,j,k)-cy))*(z(i,j,k)-cz);
//                                }
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                        {
//                            //X r = sqrt((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz));
//                            (*this)(i,j,k) = xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*0.25*sqrt(105/(2*m_pi))*((x(i,j,k)-cx)*(x(i,j,k)-cx)-(y(i,j,k)-cy)*(y(i,j,k)-cy))*(z(i,j,k)-cz);
//                        }
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalFY4(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//            X m_pi = 3.14159265358979323846;

//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitalfy4(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                {
//                                    //X r = sqrt((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz));
//                                    (*this)(i,j,k) = xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*0.25*sqrt(105/(2*m_pi))*((x(i,j,k)-cx)*(x(i,j,k)-cx)-(y(i,j,k)-cy)*(y(i,j,k)-cy))*(z(i,j,k)-cz);
//                                }
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                        {
//                            //X r = sqrt((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz));
//                            (*this)(i,j,k) = xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*0.25*sqrt(105/(2*m_pi))*((x(i,j,k)-cx)*(x(i,j,k)-cx)-(y(i,j,k)-cy)*(y(i,j,k)-cy))*(z(i,j,k)-cz);
//                        }
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalFY5(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//            X m_pi = 3.14159265358979323846;

//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitalfy5(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                {
//                                    //X r = sqrt((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz));
//                                    (*this)(i,j,k) = xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*(-0.125)*sqrt(35/m_pi)*((x(i,j,k)-cx)*(x(i,j,k)-cx)-(y(i,j,k)-cy)*(y(i,j,k)-cy))*(x(i,j,k)-cx);
//                                }
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                        {
//                            //X r = sqrt((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz));
//                            (*this)(i,j,k) = xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*(-0.125)*sqrt(35/m_pi)*((x(i,j,k)-cx)*(x(i,j,k)-cx)-(y(i,j,k)-cy)*(y(i,j,k)-cy))*(x(i,j,k)-cx);
//                        }
//#endif
//        }
//        template<class X>
//        void D3Array<X>::CalculateOrbitalFY6(const D3Array &x, const D3Array &y, const D3Array &z, const X &xs, const X &alpha, const X &cx, const X &cy, const X &cz)
//        {
//            assert (  this->m_rows == x.m_rows && this->m_columns == x.m_columns  && this->m_files == x.m_files && this->m_rows == y.m_rows && this->m_columns == y.m_columns  && this->m_files == y.m_files && this->m_rows ==z.m_rows && this->m_columns == z.m_columns  && this->m_files == z.m_files);
//            X m_pi = 3.14159265358979323846;

//#ifdef WITH_CUDA
//                if ( typeid(m_data).name() == "Pf" )
//                {
//                    scudaorbitalfy6(x.m_data,y.m_data,z.m_data,m_data,xs,cx,cy,cz,m_rows*m_columns*m_files);

//                }
//                else
//                {
//                    for(size_t i=0;i<m_rows;++i)
//                            for(size_t j=0;j<m_columns;++j)
//                                for(size_t k=0;k<m_files;++k)
//                                {
//                                    //X r = sqrt((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz));
//                                    (*this)(i,j,k) = xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*(0.125)*sqrt(35/m_pi)*((x(i,j,k)-cx)*(x(i,j,k)-cx)-(y(i,j,k)-cy)*(y(i,j,k)-cy))*(x(i,j,k)-cx);
//                                }
//                }

//#else
//            for(size_t i=0;i<m_rows;++i)
//                    for(size_t j=0;j<m_columns;++j)
//                        for(size_t k=0;k<m_files;++k)
//                        {
//                            //X r = sqrt((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz));
//                            (*this)(i,j,k) = xs*exp(-alpha*PC::AngstromToBohr*PC::AngstromToBohr*((x(i,j,k)-cx)*(x(i,j,k)-cx) + (y(i,j,k)-cy)*(y(i,j,k)-cy) + (z(i,j,k)-cz)*(z(i,j,k)-cz)))*(0.125)*sqrt(35/m_pi)*((x(i,j,k)-cx)*(x(i,j,k)-cx)-(y(i,j,k)-cy)*(y(i,j,k)-cy))*(x(i,j,k)-cx);
//                        }
//#endif
//        }

#endif //MATHTOOLS_H
