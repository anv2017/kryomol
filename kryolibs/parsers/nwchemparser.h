/*****************************************************************************************
                            nwchemparser.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef NWCHEMPARSER_H
#define NWCHEMPARSER_H

#include "parser.h"

namespace kryomol
{

class KRYOMOLPARSERS_API NwChemParser : public Parser
{
public:
  NwChemParser(const char* file);
  NwChemParser(std::istream* stream);
  ~NwChemParser();
  bool ParseFile(std::streampos pos=0);
  void ParseCouplingConstants( std::vector<QuantumCoupling>& c);
  //bool ParseFrequencies(std::streampos pos=0);
  //std::vector<Quantum::jobheader>& GetJobs();
private: 
  bool GetOptimizedGeometry();
  bool GetSinglePointGeometry();
  bool GetGeometryBlock();
  bool GetEnergies();
  int m_nsteps;
  int m_counter;
};
}
#endif
