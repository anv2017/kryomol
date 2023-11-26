/***************************************************************************
                          pdbparser.h  -  description
                             -------------------
    copyright            : (C) 2007 by A. Navarro-Vazquez
    email                : qoajnv@usc.es
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef PDBPARSER_H
#define PDBPARSER_H

#include "parser.h"

namespace nmrdev
{


/**
@author Armando Navarro-Vázquez
*/

class NMRDEVPARSERS_API PDBParser : public nmrdev::Parser
{
public:
    PDBParser(const char* file);
    PDBParser(std::istream* stream);
    ~PDBParser();
    void ParseFile(std::streampos pos=0);
private:
  std::string ExtractAtomName(const std::string& s);

};

}
#endif
