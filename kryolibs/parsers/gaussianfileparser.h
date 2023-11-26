/*****************************************************************************************
                            gaussianfileparser.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef GaussianFileParser_H
#define GaussianFileParser_H


#include <fstream>
#include <vector>

#include "parser.h"
#include "coordinate.h"
#include "frame.h"

/**
This Gaussiam16/98/03 parser will load last optimized geometry of the last job in the file
Parser implements methods for retrieving:
-Magnetic susceptibility tensors
-Scalar Couplings
-Chemical shift tensors
-ADMP BOMD dynamic runs
*/
namespace kryomol
{
class KRYOMOLPARSERS_API GaussianFileParser : public Parser
{
public:
  GaussianFileParser(const char* file);
  GaussianFileParser(std::istream* stream);
  ~GaussianFileParser();
  enum gaussversion { UNDEFINED=0,GAUSSIAN98,GAUSSIAN03,GAUSSIAN09, GAUSSIAN16 };
  bool ParseFile(std::streampos pos=0);
  bool ParseUV(std::streampos pos=0);
  virtual bool ParseFrequencies( std::streampos pos=0);
  /** Store the position of different jobs in a gaussian output*/
  std::vector<JobHeader>& Jobs();
  void ParseChemicalShifts();
  void ParseCouplingConstants(std::vector<QuantumCoupling>& c);
  D2Array<double> ParseMagneticSusceptibility();
private:
  bool ParseOrbitals(std::streampos pos);
  bool HasKeyword(std::string& line);
  bool GetGeometry();
  bool ParseArquive(std::streampos pos=0);
  bool GetVectors();
  void GetTransitionVectors(std::streampos pos);
  JobType GetJobFromRoute(const std::string& route);
  bool GetFrequencies();
  bool ExtractJBlock(std::vector<QuantumCoupling>& c);
  void GetForces(Molecule& molecule);
  bool GetDipole(const std::streampos& begin,const std::streampos& end);
  bool GetESPCharges(const std::streampos& begin, const std::streampos& end);
  bool ExistOrbitals();
  bool ExistAlphaBetaOrbitals();
  bool GetOrbitalData();
  bool GetAlphaBetaOrbitalData();
  bool GetBasisCenters();
  bool GetCoordinatesType();
  bool GetHomoLumo();
  void GetTransitionChanges();

  gaussversion GetVersion();
  kryomol::QuantumLevel GetLevel();

private:
  bool m_beta;
  int scounter;
  int m_typeD;
  int m_typeF;
  int m_homo;
  int m_lumo;
  int m_norbitals;
  gaussversion m_version;
  std::vector<std::streampos> m_pos;
  std::vector<size_t> m_atoms;
  std::vector< std::vector<Orbital> > m_orbitals;

};

}

#endif
