/*****************************************************************************************
                            hyperchemparser.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef HYPERCHEMPARSER_H
#define HYPERCHEMPARSER_H

#include "parser.h"

namespace kryomol
{
/**
A class for parsering of HyperChem files
*/
class KRYOMOLPARSERS_API HyperChemParser : public kryomol::Parser
{
public:
    HyperChemParser( const char* file);
    HyperChemParser(std::istream* stream);
    ~HyperChemParser();
    bool ParseFile(std::streampos pos=0);
private:
  void GetMolecule();
  bool GetConformer();
  bool GetConformerMolecule();
  bool m_conformation;
};

}

#endif // HYPERCHEMPARSER_H
