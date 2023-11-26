/*****************************************************************************************
                            pdbtools.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "pdbtools.h"
#include <string>

using namespace kryomol;

PDBResidue::PDBResidue(const std::string& name,const std::string& index) : _name(name), _index(index), _visible(true)
{
}

namespace kryomol
{
bool operator ==  ( const PDBResidue& a,const PDBResidue& b ) { return ( ( a.Name() == b.Name() ) && ( a.Index() == b.Index() ) ); }
}
