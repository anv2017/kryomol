/*****************************************************************************************
                            energy.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "energy.h"

using namespace kryomol;

Energy::unity Energy::m_unity=Energy::KCAL;

Energy::unity Energy::Unity() { return m_unity; }

void Energy::SetUnity(Energy::unity u) { m_unity=u; }
