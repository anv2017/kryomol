/*****************************************************************************************
                            archiveparser.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/


#include <sstream>

#include "archiveparser.h"
#include "stringtools.h"
#include "molecule.h"
#include "mathtools.h"
#include "exception.h"

using namespace kryomol;

ArchiveParser::ArchiveParser(const char* file)
    : Parser(file)
{}

ArchiveParser::ArchiveParser(std::istream* stream) : Parser(stream)
{}


ArchiveParser::~ArchiveParser()
{}

bool ArchiveParser::ParseFile(std::streampos pos)
{
  std::string archive;
  std::string line;
  m_file->clear();
  m_file->seekg(0,std::ios::beg);
  
  while(std::getline(*m_file,line) )
  {
    archive+=line;
  }
  //OK lets clear bakspaces
  std::string::size_type sit=0;
  Molecules()->push_back(Molecule());
  Molecules()->back().Frames().push_back(Frame(&Molecules()->back()));
  while( sit != std::string::npos )
  {
    sit=archive.find_first_of(" \n\t\r",sit);
    if(sit != std::string::npos)
      archive.erase(sit,1);
  }



  TokByString t(archive,"\\\\");
  std::string& geometry= t[3];

  if( t.size() < 4 ) 
  {
    throw kryomol::Exception("Error parsering gaussian archive");
  }

  StringTokenizer g1(geometry,"\\");

  Molecules()->back().Atoms().clear();
  Molecules()->back().Atoms().resize(g1.size()-1);

  size_t i;
  StringTokenizer::iterator it=g1.begin();
  it++;
  std::vector<Atom>& atoms=Molecules()->back().Atoms();
  std::vector<Coordinate>& c=Molecules()->back().Frames().back().XYZ();
  c.resize(g1.size()-1);
  for(i=0;it!=g1.end();it++,i++)
  {
    StringTokenizer g2(*it,",");
    atoms.at(i).SetSymbol(g2[0]);
    StringTokenizer::reverse_iterator it2=g2.rbegin();
    c.at(i).z()= std::atof((it2++)->c_str());
    c.at(i).y()= std::atof((it2++)->c_str());
    c.at(i).x()= std::atof(it2->c_str());
  }

  return true;

}
