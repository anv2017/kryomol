/*****************************************************************************************
                            mdlv2000parser.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef MDLV2000PARSER_H
#define MDLV2000PARSER_H

#include "parser.h"

namespace kryomol
{

class KRYOMOLPARSERS_API  MdlV2000Parser : public kryomol::Parser
{
public:
  MdlV2000Parser(const char* file);
  MdlV2000Parser(std::istream* stream);
  ~MdlV2000Parser(void);
  bool ParseFile(std::streampos pos=0);
};

}
#endif


