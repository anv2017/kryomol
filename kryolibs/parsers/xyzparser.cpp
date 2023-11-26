/*****************************************************************************************
                            xyzparser.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/


#include <iostream>
#include <sstream>
#include "xyzparser.h"
#include "molecule.h"
#include "stringtools.h"

using namespace kryomol;

XYZParser::XYZParser(const char* file)
    : Parser(file)
{}

XYZParser::XYZParser(std::istream* stream) : Parser(stream)
{}


XYZParser::~XYZParser()
{}


bool XYZParser::ParseFile(std::streampos pos)
{
  m_file->clear();
  m_file->seekg(pos,std::ios::beg);
  GetGeometry();

  return true;

}

bool XYZParser::GetGeometry()
{
  std::string line;
  Molecules()->push_back(Molecule());
  m_bfirst=false;
  while(std::getline(*m_file,line))
  {

    StringTokenizer token(line," \t\r");
    if ( token.size() == 4  )
    {
      if ( kryomol::isnum(token.at(1)) && kryomol::isnum(token.at(2)) && kryomol::isnum(token.at(3)) )
        GetMolecule(line);
    }

  }
#ifdef __GNUC__
#warning supressed move to centroid
#endif

  return true;
}

void XYZParser::GetMolecule(std::string& line)
{
  Molecule& molecule=Molecules()->back();
  molecule.Frames().push_back(Frame(&molecule));
  Frame& frame= molecule.Frames().back();
  do
  {
    StringTokenizer token(line," \t\r");
    if ( token.size() != 4 )
    {
      m_bfirst=true;
      return;
    }
    if(!m_bfirst)
    {
      Atom atom(token.at(0));
      molecule.Atoms().push_back(atom);
    }
    Coordinate c;
    c.x()=std::atof(token.at(1).c_str());
    c.y()=std::atof(token.at(2).c_str());
    c.z()=std::atof(token.at(3).c_str());

    frame.XYZ().push_back(c);

  }
  while(std::getline(*m_file,line));
}

