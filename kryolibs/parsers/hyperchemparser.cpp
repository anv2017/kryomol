/*****************************************************************************************
                            hyperchemparser.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/


#include <iostream>
#include <sstream>
#include "hyperchemparser.h"
#include "molecule.h"
#include "stringtools.h"

using namespace kryomol;

HyperChemParser::HyperChemParser(const char* file) : Parser(file)
{}

HyperChemParser::HyperChemParser(std::istream* stream) : Parser(stream)
{}


HyperChemParser::~HyperChemParser()
{}


bool HyperChemParser::ParseFile(std::streampos pos)
{
  std::string line;
  m_file->clear();
  m_file->seekg(pos,std::ios::beg);
  m_conformation = false;
  while ( std::getline ( *m_file,line ) )
  {
    if ( line.find ( "[Conformation" ) != std::string::npos )
    {
      m_conformation = true;
    }
  }
  m_file->clear();
  m_file->seekg(pos,std::ios::beg);
  while ( std::getline ( *m_file,line ) )
  {
    if ( line.find ( "mol " ) != std::string::npos )
    {
      Molecules()->push_back(Molecule());
      if (m_conformation)
          GetMolecule();
      else
          GetConformerMolecule();
    }
  }
  return true;

}

void HyperChemParser::GetMolecule()
{
  std::string line;
  Molecule& molecule=Molecules()->back();
  while ((std::getline(*m_file,line)) && (line.find ("endmol ") == std::string::npos))
  {
    if (line.find ("atom ") != std::string::npos)
    {
        StringTokenizer token(line," \t\r");
        if (token.size() >= 3)
        {
            Atom atom(token.at(3));
            molecule.Atoms().push_back(atom);
        }
    }
  }
  while ( std::getline(*m_file,line) /* && (line.find ("[HIN System Description]) == std::string::npos)*/)
  {
    if (line.find ("[Conformation") != std::string::npos)
    {
        GetConformer();
    }
  }
}

bool HyperChemParser::GetConformer()
{
  std::string line;
  Molecule& molecule=Molecules()->back();
  molecule.Frames().push_back(Frame(&molecule));
  Frame& frame= molecule.Frames().back();
  while ((std::getline(*m_file,line)) && (line.find ("X(") == std::string::npos))
  {
    if (line.find("Energy") != std::string::npos)
    {
        StringTokenizer token ( line,"=" );
        if ( kryomol::isnum(token.at(1)))
        {
            frame.PotentialEnergy() = kryomol::atof(token.at(1));
        }
    }
  }
  do
  {
    StringTokenizer token(line,"=");
    line = token.at(1);
    token=StringTokenizer(line," \t\r");
    if ( kryomol::isnum(token.at(0)) && kryomol::isnum(token.at(1)) && kryomol::isnum(token.at(2)) )
    {
        Coordinate c;
        c.x()=kryomol::atof(token.at(0));
        c.y()=kryomol::atof(token.at(1));
        c.z()=kryomol::atof(token.at(2));

        frame.XYZ().push_back(c);
    }
  }
  while (std::getline(*m_file,line) && (line.find ("X(") != std::string::npos));

#ifdef __GNUC__
#warning supressed move to centroid
#endif

  return true;
}

bool HyperChemParser::GetConformerMolecule()
{
  std::string line;
  Molecule& molecule=Molecules()->back();
  molecule.Frames().push_back(Frame(&molecule));
  Frame& frame= molecule.Frames().back();
  while ((std::getline(*m_file,line)) && (line.find ("endmol ") == std::string::npos))
  {
    if (line.find ("atom ") != std::string::npos)
    {
        StringTokenizer token(line," \t\r");
        if (token.size() >= 10)
        {
            Atom atom(token.at(3));
            molecule.Atoms().push_back(atom);
            if ( kryomol::isnum(token.at(7)) && kryomol::isnum(token.at(8)) && kryomol::isnum(token.at(9)) )
            {
                Coordinate c;
                c.x()=kryomol::atof(token.at(7));
                c.y()=kryomol::atof(token.at(8));
                c.z()=kryomol::atof(token.at(9));

                frame.XYZ().push_back(c);
            }
        }
    }
  }
#ifdef __GNUC__
#warning supressed move to centroid
#endif

  return true;
}


