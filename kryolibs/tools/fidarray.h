/*****************************************************************************************
                            fidarray.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef FIDARRAY_H
#define FIDARRAY_H

#include <valarray>
#include <complex>


/**
  * container adapted to do basic NMR processing operations
  * FFT code adapted from the GNU scientific library (gsl)
  * as from its version 1.3
  */

template <class T> inline std::complex<T> expi(float a)
{
  return std::complex<T>(cos(a),sin(a));
}

template <class T> inline std::complex<T> expi(double a)
{
  return std::complex<T>(cos(a),sin(a));
}

  /** implemented by similarity to vector */
class fidarray : public std::valarray< std::complex<float> >
{
public:
  fidarray();
  fidarray(size_t n);
  fidarray& shift(int np);
  fidarray& cshift(int np);
  ~fidarray();

  /** Same behaviour as resize in std::vector */
  void grow(size_t n);
  /** is array empty*/
  bool empty() const{ return (this->size() == 0 );}
  /** multiply odd points by -1 */
  fidarray& quadrature();

  /** implemented by similarity to vector */
  void clear();

  /** fast fourier transform*/
  void fft(int sense);
  /** wavelet wlt transform*/
  void wlt(int sense, int coefficients);
  /** real FFT */
  void fft_real();
  /** returns the total sum of real or imaginary points*/
  float sum(bool real=true);
  /** sum of squares */
  void magnitude();

  /** squared root of sum of squares */
  void power();

  /** Do Bit reversion for Fourier transform,
    * adapted from the GNU scientific library*/
  void bitreverse();

  /** Do Bit reversion for real Fourier transform,
    * adapted from the GNU scientific library*/
  void bitreverse_real();

  /** Get the log2(size) for fft
    * adapted from the GNU scientific library*/
  size_t logsize(bool& ispoweroftwo);
  /** convert a real (sequential bruker) spectrum into a normal complex one */
  void real2complex();
  /** divide all values by the size*/
  void normalize();

};

#endif
