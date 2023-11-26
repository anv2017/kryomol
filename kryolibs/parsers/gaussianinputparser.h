/*****************************************************************************************
                            gaussianinputparser.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef GAUSSIANINPUTPARSER_H
#define GAUSSIANINPUTPARSER_H

#include "parser.h"

namespace kryomol
{
/**
A Gaussian98/03 inoput file parser
*/
class KRYOMOLPARSERS_API GaussianInputParser : public kryomol::Parser
{
public:
  GaussianInputParser(const char* file);
  GaussianInputParser(std::istream* stream);
  ~GaussianInputParser();
  bool ParseFile(std::streampos pos=0);
private:
  void GetGeometry();

};
}
#endif
