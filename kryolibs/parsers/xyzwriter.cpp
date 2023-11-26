/*****************************************************************************************
                            xyzwriter.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <iostream>

#include  "molecule.h"
#include "xyzwriter.h"

using namespace kryomol;

XYZWriter::XYZWriter(Molecule* molecule, bool all) : m_molecule(molecule), m_all(all)
{
}

XYZWriter::~XYZWriter()
{
}

std::ostream& kryomol::operator << (std::ostream& s, const XYZWriter& w)
{
    std::cout << "Exporting in XYZ format" << std::endl;
    s << w.m_molecule->Atoms().size() << std::endl;
    if (w.m_all)
    {
        for (size_t i=0; i<w.m_molecule->Frames().size();++i)
        {
            s << "Frame #" << i+1 << std::endl;
            s << w.m_molecule->Frames().at(i) << std::endl;
            s << std::endl;
        }

    }
    else
    {
        s << "title" << std::endl;
        s << w.m_molecule->CurrentFrame() << std::endl;
    }
    return s;

}


