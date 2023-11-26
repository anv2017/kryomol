/*****************************************************************************************
                            cpmdwriter.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef CPMDWRITER_H
#define CPMDWRITER_H

#include <iostream>
#include "molecule.h"
#include "stringtools.h"
#include "parsersexport.h"
#include <sstream>
#include <vector>
#include <iomanip>

/**
A class to write basic cpmd input files
*/

namespace kryomol{

    class Molecule;

class KRYOMOLPARSERS_API CPMDWriter
{
public:
    CPMDWriter(kryomol::Molecule* molecule);
    ~CPMDWriter();

    friend KRYOMOLPARSERS_API void operator << (std::ostream& s, const CPMDWriter& aw);

private:
    kryomol::Molecule* m_molecule;

};

}

#endif
