/*****************************************************************************************
                            archiveparser.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef ARCHIVEPARSER_H
#define ARCHIVEPARSER_H

#include "parser.h"

namespace kryomol
{
/**
Gaussian Archive Parser.  It can extract geometry forces and hessian
*/

class KRYOMOLPARSERS_API ArchiveParser : public kryomol::Parser
{
public:
    ArchiveParser(const char* file);
    ArchiveParser(std::istream* stream);
    ~ArchiveParser();    
    bool ParseFile(std::streampos pos=0 );
};

}
#endif
