/*****************************************************************************************
                            gamesparser.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef GAMESSPARSER_H
#define GAMESSPARSER_H

#include "parser.h"

namespace kryomol
{

class KRYOMOLPARSERS_API GamessParser : public Parser
{
public:
    GamessParser( const char* file);
    GamessParser(std::istream* stream);
    ~GamessParser();
    bool ParseFile(std::streampos pos=0);
    void ParseChemicalShifts();
    //virtual std::vector<Quantum::jobheader>& GetJobs();

private:
  void GetGeometry();
 bool ExtractTensorBlock();
};

}
#endif
