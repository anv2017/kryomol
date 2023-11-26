/*****************************************************************************************
                            orbital.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef ORBITAL_H
#define ORBITAL_H

#include <iostream>
#include <vector>
#include "mathtools.h"
#include "coreexport.h"
#include "coordinate.h"

namespace kryomol
{
    class Coordinate;

/** @brief atomic and molecular orbitals data

This class manage atomic and molecular orbitals data from compiled from parser file*/
class KRYOMOLCORE_API Orbital
{
public:
    enum OrbitalType {S, SP, P, D, F};

    Orbital();
    Orbital(OrbitalType orbitaltype, std::vector<float> alpha, std::vector<float> xs, std::vector<float> xp) { m_orbitaltype=orbitaltype;  m_alpha=alpha; m_xs=xs; m_xp=xp;}

    OrbitalType Type() { return m_orbitaltype; }
    std::vector<float>& Xs() { return m_xs; }
    std::vector<float>& Xp() { return m_xp; }
    std::vector<float>& Alpha() { return m_alpha; }

private:
    Coordinate m_atom;
    OrbitalType m_orbitaltype;
    std::vector<float> m_xs;
    std::vector<float> m_xp;
    std::vector<float> m_alpha;

};

}

#endif // ORBITAL_H
