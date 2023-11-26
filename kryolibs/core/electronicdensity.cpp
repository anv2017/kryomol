/*****************************************************************************************
                            electronicdensity.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "electronicdensity.h"

using namespace kryomol;

ElectronicDensity::ElectronicDensity()
{}

ElectronicDensity::ElectronicDensity(int nx, int ny, int nz, float dx, float dy, float dz, const Coordinate &origin) : m_nx(nx), m_ny(ny), m_nz(nz), m_dx(dx), m_dy(dy), m_dz(dz), m_origin(origin)
{}

