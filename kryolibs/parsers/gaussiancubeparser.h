/*****************************************************************************************
                            gaussiancubeparser.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef GAUSSIANCUBEPARSER_H
#define GAUSSIANCUBEPARSER_H

#include "parser.h"

namespace kryomol
{
/**
A class for parsering of GaussianCube files
*/
class KRYOMOLPARSERS_API GaussianCubeParser : public kryomol::Parser
{
public:
    GaussianCubeParser( const char* file);
    GaussianCubeParser(std::istream* stream);
    ~GaussianCubeParser();
    bool ParseFile(std::streampos pos=0);
private:
};

}

#endif // GAUSSIANCUBEPARSER_H
