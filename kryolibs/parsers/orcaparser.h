/*****************************************************************************************
                            orcaparser.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef ORCAPARSER_H
#define ORCAPARSER_H

#include <fstream>
#include <vector>

#include "parser.h"
#include "coordinate.h"
#include "frame.h"
#include "stringtools.h"


namespace kryomol
{
class KRYOMOLPARSERS_API OrcaParser : public Parser
{
public:
  OrcaParser(const char* file);
  OrcaParser(std::istream* stream);
  ~OrcaParser();
  bool ParseFile(std::streampos pos=0);
  std::vector<JobHeader>& Jobs();
  bool ParseUV(std::streampos pos=0);
  bool ParseFrequencies(std::streampos pos=0);
private:
  bool ParseOrbitals(std::streampos pos);
  bool GetGeometry();
  bool ExistOrbitals();
  bool GetOrbitalData();
  bool GetBasisCenters();
  bool GetHomoLumo();
  void GetEnergyForFrame();
  void GetGradientForFrame();
  void ParseUVLengthBlock(std::vector<Spectralline>& lines);
  void ParseUVVelocityBlock(std::vector<Spectralline>& lines);
  void ParseCDLengthBlock(std::vector<Spectralline>& lines);
  void ParseCDVelocityBlock(std::vector<Spectralline>& lines);
  void ParseSolventShiftBlock(std::vector<Spectralline>& lines);
  bool GetFrequencies();
  bool GetNormalModes();
  bool GetNormalModeBlock(StringTokenizer& headertok);


private:
  int m_homo;
  int m_lumo;
  int m_norbitals;
  std::vector<std::streampos> m_pos;
  std::vector<size_t> m_atoms;
  std::vector< std::vector<Orbital> > m_orbitals;
  std::vector < std::string > m_jobkeys;

};

}

#endif // ORCAPARSER_H
