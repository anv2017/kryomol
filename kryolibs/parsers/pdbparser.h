/*****************************************************************************************
                            pdbparser.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef PDBPARSER_H
#define PDBPARSER_H

#include <vector>
#include <string>

#include "parser.h"

namespace kryomol
{

class KRYOMOLPARSERS_API PDBParser : public kryomol::Parser
{
public:
    PDBParser(const char* file);
    PDBParser(std::istream* stream);
    ~PDBParser();
    bool ParseFile(std::streampos pos=0);
private:
  std::vector<std::string> ParseAtomLine(const std::string& line);
  std::string ExtractAtomName(const std::string& s);

};

}
#endif
