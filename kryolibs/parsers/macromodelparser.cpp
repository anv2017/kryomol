/*****************************************************************************************
                            macromodelparser.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/


#include <sstream>
#include "macromodelparser.h"
#include "stringtools.h"
#include "molecule.h"
#include "exception.h"

using namespace kryomol;

MacroModelParser::MacroModelParser(const char* file) : Parser(file)
{}

MacroModelParser::MacroModelParser(std::istream* stream) : Parser(stream)
{}

bool MacroModelParser::ParseFile(std::streampos pos)
{
  m_file->clear();
  m_file->seekg(0,std::ios::beg);
  std::string line;
  std::getline(*m_file,line);
  //Get the number of atoms
  StringTokenizer tok(line," \t=");
  int natoms=atoi(tok.front().c_str());
  StringTokenizer::iterator tt;
  double venergy=0;
  for(tt=tok.begin();tt!=tok.end();++tt)
  {
    if ( *tt == "E" ) venergy=std::atof((tt+1)->c_str());
  }


  if( natoms < 1) 
  {
    throw kryomol::Exception("Error parsering macromodel file");
  }

  Molecules()->push_back(Molecule());
  Molecule& molecule=Molecules()->back();
  molecule.Atoms().resize(natoms);
  molecule.Frames().push_back(Frame(&Molecules()->back()));
  Frame& frame=molecule.Frames().back();
  frame.XYZ().resize(natoms);

  size_t i=0;
  //Setup first structure
  while(std::getline(*m_file,line))
  {
    StringTokenizer token(line," \t");
    if(token.size() < 10 ) break;
    molecule.Atoms().at(i).SetSymbol(GetAtom(atoi(token.front().c_str())));
    std::stringstream str;
    str << i+1;
    frame.XYZ().at(i).x()=atof(token.at(13));
    frame.XYZ().at(i).y()=atof(token.at(14));
    frame.XYZ().at(i).z()=atof(token.at(15));
    i++;
  }
  frame.KineticEnergy()=0;
  frame.PotentialEnergy()=0.2390*venergy;

  bool hasconexions=true;
  do
  {
    molecule.Frames().push_back(Frame(&molecule));
    Frame& frame=molecule.Frames().back();
    frame.XYZ().resize(natoms);
    
    //Get Energies
    StringTokenizer tok(line," \t=");
    if (tok.empty() )  break;
    StringTokenizer::iterator tt;
    if ( atoi(tok.front().c_str() ) < 0 ) hasconexions=false;
    for(tt=tok.begin();tt!=tok.end();++tt)
    {
      if ( *tt == "KE" )
        frame.KineticEnergy()=0.2390*std::atof((tt+1)->c_str());
      if ( *tt == "PE" || *tt == "E" )
        frame.PotentialEnergy()=0.2390*std::atof((tt+1)->c_str());
    }
    i=0;
    while(std::getline(*m_file,line))
    {

      StringTokenizer tok(line," \t");

      if(!hasconexions)
      {
        frame.XYZ().at(i).x()=std::atof(tok.at(1).c_str());
        frame.XYZ().at(i).y()=std::atof(tok.at(2).c_str());
        frame.XYZ().at(i).z()=std::atof(tok.at(3).c_str());
      }
      else
      {
        frame.XYZ().at(i).x()=std::atof(tok.at(13).c_str());
        frame.XYZ().at(i).y()=std::atof(tok.at(14).c_str());
        frame.XYZ().at(i).z()=std::atof(tok.at(15).c_str());
      }
      i++;
      if( i == (size_t) natoms ) break;
    }

  }
  while (std::getline(*m_file,line));

  return true;


}


MacroModelParser::~MacroModelParser()
{}

int MacroModelParser::GetBondType(int bt)
{
  switch(bt)
  {
  case 0:
    return 0;
  case 1:
  case 4:
    return 1;
  case 2:
    return 2;
  case 3:
    return 3;
  default:
    return 0;

  }
}
std::string MacroModelParser::GetAtom(int type)
{
  switch(type)
  {
  case 1:
  case 2:
  case 3:
  case 14:
    return "C";
  case 15:
  case 16:
  case 18:
  case 23:
    return "O";
  case 24:
  case 25:
  case 26:
  case 31:
  case 32:
  case 40:
    return "N";
  case 41:
  case 42:
  case 43:
  case 44:
  case 45:
  case 48:
    return "H";
  case 49:
  case 51:
  case 52:
    return "S";
  case 53:
    return "P";
  case 56:
    return "F";
  case 57:
    return "Cl";
  case 58:
    return "Br";
  case 59:
    return "I";
  case 60:
    return "Si";
  default:
    return "X";

  }
}


