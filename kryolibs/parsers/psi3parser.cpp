/***************************************************************************
 *   Copyright (C) 2004 by Armando Navarro Vázquez                         *
 *   armando@rmntec                                                        *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <iostream>
#include <sstream>


#include "psi3parser.h"
#include "tokenizer.h"
#include "molecule.h"

using namespace nmrdev;
const double bohrtoangs=0.529189379;

Psi3Parser::Psi3Parser(const char* file)
 : Parser(file)
{
}

Psi3Parser::Psi3Parser(std::istream* stream) : Parser(stream)
{}

Psi3Parser::~Psi3Parser()
{
}

bool Psi3Parser::ParseFile(std::streampos pos/*0*/)
{
  std::cout << "Psi3Parser" << std::endl;
  std::string	line;
  
  std::vector<std::streampos> steps;
  while(std::getline(*m_file,line))
  {
     std::streampos pos=m_file->tellg();
     if(line.find("Cartesian geometry and possibly gradient",0)!=std::string::npos )
     {
       steps.push_back(pos);
     }
  }
  m_file->clear();
  m_file->seekg(0,std::ios::beg);
  std::vector<std::streampos>::iterator it;
  for(it=steps.begin();it!=steps.end();it++)
  {
    CMolecule molecule;
    m_file->seekg(*it,std::ios::beg);
    while(std::getline(*m_file,line))
    {
    
     StringTokenizer token(line," \t");
     size_t size=token.size();
     std::stringstream label;
     int counter=0;
     if(token.size() == 5)
     {
       counter++;
       CAtom atom;
       atom.c.x=bohrtoangs*std::atof(token[size-3].c_str());
       atom.c.y=bohrtoangs*std::atof(token[size-2].c_str());
       atom.c.z=bohrtoangs*std::atof(token[size-1].c_str());
       label << counter;
       atom.label=label.str();
       atom.name=CMolecule::GetAtomName(std::atoi(token[0].c_str()));
       molecule.m_atom.push_back(atom);
     }
     else
       break;
    }
    molecule.SetConexions();
    Molecules()->push_back(molecule);
  }

  return true;
}


