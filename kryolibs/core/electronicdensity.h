/*****************************************************************************************
                            electronicdensity.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef ELECTRONICDENSITY_H
#define ELECTRONICDENSITY_H

#include <iostream>

#include "mathtools.h"
#include "orbitalarray.h"
#include "coreexport.h"
#include "coordinate.h"

namespace kryomol
{

/** @brief electronic density

This class manage the electronic density of a molecule*/
class KRYOMOLCORE_API ElectronicDensity
{
public:
    enum unity { ANGSTROM, BOHR };

    ElectronicDensity();
    ElectronicDensity(int nx, int ny, int nz, float dx, float dy, float dz, const Coordinate& origin=Coordinate(0,0,0));
    ~ElectronicDensity() {}

    int Nx() const { return m_nx; }
    int Ny() const { return m_ny; }
    int Nz() const { return m_nz; }
    float Dx() const { return m_dx; }
    float Dy() const { return m_dy; }
    float Dz() const { return m_dz; }
    const OrbitalArray& Density() const {return m_density;}
    OrbitalArray& Density() {return m_density;}
    const Coordinate& Origin() const {return m_origin;}
    Coordinate& Origin() {return m_origin;}

    void SetDensity (const OrbitalArray& density) {m_density = density;}
    void SetOrigin(const Coordinate& origin) {m_origin = origin;}

private:
    int m_nx;
    int m_ny;
    int m_nz;
    float m_dx;
    float m_dy;
    float m_dz;
    Coordinate m_origin;
    OrbitalArray m_density;

};

}

#endif // ELECTRONICDENSITY_H
