/*****************************************************************************************
                            macromodelparser.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef MACROMODELPARSER_H
#define MACROMODELPARSER_H
#include "parser.h"


namespace kryomol
{
/**
A parser for macromodel molecular dynamics
*/

class KRYOMOLPARSERS_API MacroModelParser : public kryomol::Parser
{
public:
    MacroModelParser(const char*file);
    MacroModelParser(std::istream* stream);
    ~MacroModelParser();
    bool ParseFile(std::streampos pos=0);
private:
    std::string GetAtom(int type);
    int GetBondType(int bt);

};

}
#endif
