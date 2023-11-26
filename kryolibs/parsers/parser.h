/*****************************************************************************************
                            parser.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QUANTUMPARSER_H
#define QUANTUMPARSER_H

#include <fstream>
#include <istream>
#include <vector>
#include "parsersexport.h"

#include "quantumcoupling.h"


namespace kryomol
{
  class Molecule;
  enum QuantumLevel { MNDO, HF, MP2, QCISD, CCSD, CCSDT, SCRF };
  enum JobType {singlepoint, opt, freq, uv, nmr, dyn};

  /** A simple structure to store the position of different jobs in gaussian files*/
  struct JobHeader
  {
    JobHeader ( JobType t, std::streampos p ) : type ( t ) , pos ( p ) {}
    JobType type;
    std::streampos pos;
  };

  /**
  The base class for all molecular structure parsers
  */

  class KRYOMOLPARSERS_API Parser
  {
    public:
      Parser ( const char* inputfile );
      Parser ( std::istream* stream );
      virtual ~Parser();
      void SetMolecules ( std::vector<kryomol::Molecule>* molecule ) { m_molecules=molecule; }
      void Parse(std::streampos pos=0);     
      virtual bool ParseFile(std::streampos pos=0 )=0;
      //virtual void ParseFile ( std::streampos pos) =0;
      virtual void ParseChemicalShifts();
      virtual void ParseCouplingConstants(std::vector<QuantumCoupling>& c);
      virtual D2Array<double> ParseMagneticSusceptibility();
      virtual bool ParseFrequencies( std::streampos =0) { return false; }
      virtual bool ParseUV( std::streampos =0) { return false; }
      virtual std::vector<JobHeader>& Jobs();
    protected:
      const std::vector<kryomol::Molecule>* Molecules() const { return m_molecules; }
      std::vector<kryomol::Molecule>* Molecules() { return m_molecules; }
    protected:
      std::istream* m_file;
      QuantumLevel m_level;
      std::vector<JobHeader> m_jobpos;
    private:
      std::vector<kryomol::Molecule>* m_molecules;
      bool m_bcreated;

  };

};
#endif
