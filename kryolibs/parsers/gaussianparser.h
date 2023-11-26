/*****************************************************************************************
                            gaussianparser.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef GAUSSIANPARSER_H
#define GAUSSIANPARSER_H


#include <fstream>
#include <vector>

#include "parser.h"

/**
This Gaussiah98/03 parser will load last optimized geometry of the last job in the file
Parser implements methods for retrieving:
-Magnetic susceptibility tensors
-Scalar Couplings
-Chemical shift tensors
*/
namespace kryomol
{
class KRYOMOLPARSERS_API GaussianParser : public Parser
{
public:
  GaussianParser(const char* file);
  GaussianParser(std::istream* stream);
  ~GaussianParser();
  bool ParseFile(std::streampos pos=0);
  /** Store the position of different jobs in a gaussian output*/
  std::vector<JobHeader>& Jobs();
  void ParseChemicalShifts();
  void ParseCouplingConstants(std::vector<QuantumCoupling>& c);
  D2Array<double> ParseMagneticSusceptibility();
private:
  bool HasKeyword(std::string& line);
  bool GetGeometry();
  bool ExtractJBlock(std::vector<QuantumCoupling>& c);
  int scounter;

};

}

#endif
