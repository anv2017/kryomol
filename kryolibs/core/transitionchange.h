/*****************************************************************************************
                            orbitaldata.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef TRANSITIONCHANGE_H
#define TRANSITIONCHANGE_H

#include <iostream>
#include <vector>
#include "mathtools.h"
#include "coreexport.h"

namespace kryomol
{

/** @brief transition change data

This class manage transition change data */
class KRYOMOLCORE_API TransitionChange
{
public:
    TransitionChange();
    TransitionChange(int i, int j, float coeff);
    TransitionChange(std::string& i, std::string& j, float coeff);
    ~TransitionChange() {}

    int OrbitalI() const { return m_i;}
    int OrbitalJ() const {return m_j; }
    std::string& OrbitalSI() {return m_si;}
    const std::string& OrbitalSI() const {return m_si;}
    std::string& OrbitalSJ() {return m_sj;}
    const std::string& OrbitalSJ() const {return m_sj;}
    float Coefficient() const { return m_coefficient; }

private:
    int m_i;
    int m_j;
    std::string m_si;
    std::string m_sj;
    float m_coefficient;

};

}

#endif // TRANSITIONCHANGE_H
