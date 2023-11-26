/*****************************************************************************************
                            xyzwriter.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef XYZWRITER_H
#define XYZWRITER_H

#include "parsersexport.h"
#include "sstream"

namespace kryomol
{
  class Molecule;

class KRYOMOLPARSERS_API XYZWriter
{
public:
    XYZWriter(kryomol::Molecule* molecule, bool all);
    ~XYZWriter();
    friend KRYOMOLPARSERS_API std::ostream& operator << (std::ostream& s, const XYZWriter& aw);
private:
   kryomol::Molecule* m_molecule;
   bool m_all;

};
  
}

#endif
