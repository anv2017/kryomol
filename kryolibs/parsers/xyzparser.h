/*****************************************************************************************
                            xyzparser.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef XYZPARSER_H
#define XYZPARSER_H

#include "parser.h"

namespace kryomol
{
/**
A class for parsering of XYZ files
*/
class KRYOMOLPARSERS_API XYZParser : public kryomol::Parser
{
public:
    XYZParser( const char* file);
    XYZParser(std::istream* stream);
    ~XYZParser();
    bool ParseFile(std::streampos pos=0);
private:
  bool GetGeometry();
  void GetMolecule(std::string& line);
  bool m_bfirst;

};

}

#endif
