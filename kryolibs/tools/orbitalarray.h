/*****************************************************************************************
                            orbitalarray.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef ORBITALARRAY_H
#define ORBITALARRAY_H

#include <math.h>
#include "mathtools.h"

namespace kryomol
{

class OrbitalArray : public D3Array<float>
{
public:
    /** Build an empty orbital array */
    OrbitalArray():D3Array<float>(){}
    /** Build an orbital array with  i rows and j columns*/
    explicit OrbitalArray ( size_t i, size_t j, size_t k ):D3Array<float>( i, j, k ){}
    /** Build an orbital array with  i rows and j columns*, with values x*/
    OrbitalArray ( size_t i,size_t j, size_t k, const float x ):D3Array<float>(i,j,k,x){}
    /** copy constructor*/
    OrbitalArray ( const OrbitalArray& x ):D3Array<float>(x){}
    /** destructor, erase all elements from the heap*/
    ~OrbitalArray(){}
    /** calculate the array of exponentails used for calculating the atomic orbitals*/
    void CalculateExponential(const float Nx, const float Ny, const float Nz, const float step, const float cx, const float cy, const float cz, const float alpha, const float xs);
    /** calculate the orbital type S*/
    void CalculateOrbitalS (  const OrbitalArray& exponential);
    /** calculate the orbital type P*/
    void CalculateOrbitalPx (  const OrbitalArray& exponential, const float cx, const float Nx, const float step);
    void CalculateOrbitalPy (  const OrbitalArray& exponential, const float cy, const float Ny, const float step);
    void CalculateOrbitalPz (  const OrbitalArray& exponential, const float cz, const float Nz, const float step);
    /** calculate the orbital type D*/
    void CalculateOrbitalDxx (  const OrbitalArray& exponential, const float cx, const float Nx, const float step);
    void CalculateOrbitalDyy (  const OrbitalArray& exponential, const float cy, const float Ny, const float step);
    void CalculateOrbitalDzz (  const OrbitalArray& exponential, const float cz, const float Nz, const float step);
    void CalculateOrbitalDxy (  const OrbitalArray& exponential, const float cx, const float cy, const float Nx, const float Ny, const float step);
    void CalculateOrbitalDxz (  const OrbitalArray& exponential, const float cx, const float cz, const float Nx, const float Nz, const float step);
    void CalculateOrbitalDyz (  const OrbitalArray& exponential, const float cy, const float cz, const float Ny, const float Nz, const float step);
    /** calculate the orbital type F*/
    void CalculateOrbitalFxxx (  const OrbitalArray& exponential, const float cx, const float Nx, const float step);
    void CalculateOrbitalFyyy (  const OrbitalArray& exponential, const float cy, const float Ny, const float step);
    void CalculateOrbitalFzzz (  const OrbitalArray& exponential, const float cz, const float Nz, const float step);
    void CalculateOrbitalFxxy (  const OrbitalArray& exponential, const float cx, const float cy, const float Nx, const float Ny, const float step);
    void CalculateOrbitalFxyy (  const OrbitalArray& exponential, const float cx, const float cy, const float Nx, const float Ny, const float step);
    void CalculateOrbitalFxxz (  const OrbitalArray& exponential, const float cx, const float cz, const float Nx, const float Nz, const float step);
    void CalculateOrbitalFxzz (  const OrbitalArray& exponential, const float cx, const float cz, const float Nx, const float Nz, const float step);
    void CalculateOrbitalFyyz (  const OrbitalArray& exponential, const float cy, const float cz, const float Ny, const float Nz, const float step);
    void CalculateOrbitalFyzz (  const OrbitalArray& exponential, const float cy, const float cz, const float Ny, const float Nz, const float step);
    void CalculateOrbitalFxyz (  const OrbitalArray& exponential, const float cx, const float cy, const float cz, const float Nx, const float Ny, const float Nz, const float step);
    /** calculate the orbital type D with spherical harmonics*/
    void CalculateOrbitalDY0 (  const OrbitalArray& exponential, const float cx, const float cy, const float cz, const float Nx, const float Ny, const float Nz, const float step);
    void CalculateOrbitalDY1 (  const OrbitalArray& exponential, const float cx, const float cz, const float Nx, const float Nz, const float step);
    void CalculateOrbitalDY2 (  const OrbitalArray& exponential, const float cy, const float cz, const float Ny, const float Nz, const float step);
    void CalculateOrbitalDY3 (  const OrbitalArray& exponential, const float cx, const float cy, const float Nx, const float Ny, const float step);
    void CalculateOrbitalDY4 (  const OrbitalArray& exponential, const float cx, const float cy, const float Nx, const float Ny, const float step);
    /** calculate the orbital type F with spherical harmonics*/
    void CalculateOrbitalFY0 (  const OrbitalArray& exponential, const float cx, const float cy, const float cz, const float Nx, const float Ny, const float Nz, const float step);
    void CalculateOrbitalFY1 (  const OrbitalArray& exponential, const float cx, const float cy, const float cz, const float Nx, const float Ny, const float Nz, const float step);
    void CalculateOrbitalFY2 (  const OrbitalArray& exponential, const float cx, const float cy, const float cz, const float Nx, const float Ny, const float Nz, const float step);
    void CalculateOrbitalFY3 (  const OrbitalArray& exponential, const float cx, const float cy, const float cz, const float Nx, const float Ny, const float Nz, const float step);
    void CalculateOrbitalFY4 (  const OrbitalArray& exponential, const float cx, const float cy, const float cz, const float Nx, const float Ny, const float Nz, const float step);
    void CalculateOrbitalFY5 (  const OrbitalArray& exponential, const float cx, const float cy, const float Nx, const float Ny, const float step);
    void CalculateOrbitalFY6 (  const OrbitalArray& exponential, const float cx, const float cy, const float Nx, const float Ny, const float step);

};

}

#endif // ORBITALARRAY_H
