/*****************************************************************************************
                            transitionchange.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "transitionchange.h"

using namespace kryomol;

TransitionChange::TransitionChange()
{}

TransitionChange::TransitionChange(int i, int j, float coeff)
{
    m_i = i;
    m_j = j;
    m_coefficient = coeff;
}

TransitionChange::TransitionChange(std::string& i, std::string& j, float coeff)
{
    m_si = i;
    m_sj = j;
    m_coefficient = coeff;
}
