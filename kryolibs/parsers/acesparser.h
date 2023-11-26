/*****************************************************************************************
                            acesparser.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef ACESPARSER_H
#define ACESPARSER_H

#include "parser.h"

namespace kryomol
{
/**
@author Armando Navarro-Vazquez
*/
class KRYOMOLPARSERS_API AcesParser : public kryomol::Parser
{
public:
    AcesParser(const char* file);
    AcesParser(std::istream* stream);
    ~AcesParser();

    bool ParseFile(std::streampos pos=0);
    bool ParseFrequencies(std::streampos pos=0);
    std::vector<JobHeader>& Jobs();

private:
  bool GetGeometry();
  void GetFrequencies();
  void GetVectors();
  bool GetVectorBlock(int nmodes);
  kryomol::QuantumLevel GetLevel();
  int m_counter;

};
}

#endif
