/*****************************************************************************************
                            aceswriter.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef ACESWRITER_H
#define ACESWRITER_H

#include "stringtools.h"
#include "parsersexport.h"

#include "molecule.h"
#include <sstream>
#include <vector>
#include <iomanip>

namespace kryomol{

class Molecule;

class KRYOMOLPARSERS_API AcesWriter
{
public:
    AcesWriter(kryomol::Molecule* molecule);// : m_molecule(molecule) {}
    ~AcesWriter();

    friend KRYOMOLPARSERS_API std::ostream& operator << (std::ostream& s, const AcesWriter& aw);
    //friend KRYOPARSERS_EXPORT std::ostream& operator << (std::ostream& s, const AcesWriter& aw);


private:
    kryomol::Molecule* m_molecule;
    struct zmatlabel { std::stringstream label; float value; };

};


}


#endif
