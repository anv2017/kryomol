/*****************************************************************************************
                            basiscenter.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef BASISCENTER_H
#define BASISCENTER_H

#include <vector>
#include "orbital.h"

namespace kryomol {

    class Coordinate;

/** @brief basis functions for atomic and molecular orbitals

This class manage the basis functions for making atomic and molecular orbitals*/

class KRYOMOLCORE_API BasisCenter
{
public:
    BasisCenter();
    BasisCenter(Coordinate& atom, std::vector<Orbital>& orbitals) {m_atom = atom; m_orbitals=orbitals;}
    const Coordinate& Atom() const {return m_atom;}
    Coordinate& Atom() {return m_atom;}
    const std::vector<Orbital>& Orbitals() const {return m_orbitals;}
    std::vector<Orbital>& Orbitals() {return m_orbitals;}

private:
    std::vector<Orbital> m_orbitals;
    Coordinate m_atom;
};

}

#endif // BASISCENTER_H
