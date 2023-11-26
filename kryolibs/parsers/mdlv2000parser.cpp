/*****************************************************************************************
                            mdlv2000parser.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "mdlv2000parser.h"
#include "molecule.h"
#include "stringtools.h"
#include "exception.h"
#include <vector>
#include <sstream>

using namespace kryomol;

MdlV2000Parser::MdlV2000Parser ( const char* file ) : Parser ( file )
{}

MdlV2000Parser::MdlV2000Parser ( std::istream* stream ) : Parser ( stream )
{}

MdlV2000Parser::~MdlV2000Parser ( void )
{}


bool MdlV2000Parser::ParseFile ( std::streampos pos )
{
  std::string line;
  m_file->clear();
  m_file->seekg ( pos,std::ios::beg );
  bool bfirst=true;
  Molecules()->push_back ( Molecule() );
  while ( std::getline ( *m_file,line ) )
  {

    //skip three first lines
    size_t i=0;


    //OK, gets the number of atoms and bonds

    while ( std::getline ( *m_file,line ) )
    {
      StringTokenizer token ( line," \t\r" );
      size_t size=token.size();
      if ( size )
        if ( size == 11L || token.back() == "V2000" ) goto lb1;
    }
    return true;
  lb1:
    StringTokenizer token ( line," \t\r" );

    if ( token.empty() )
    {
      throw kryomol::Exception ( "Error parsering mol file" );
    }
    //either there are 11 numbers or last token is "V2000"
    if ( token.size() != 11L  && token.back() !="V2000" )
    {
      throw kryomol::Exception ( "Error parsering mdl V2000 file" );
    }
    bfirst=false;
    int natoms=atoi ( token.at(0).c_str() );
    if ( natoms < 1 )
    {
      throw kryomol::Exception ( "Error parsering mdl V2000 file" );
    }
    Frame frame ( &Molecules()->back() );
    Molecules()->back().Frames().push_back ( frame );

    Molecules()->back().Atoms().clear();
    Molecules()->back().Atoms().resize ( natoms );
    Molecules()->back().Frames().back().XYZ().resize ( natoms );

    std::vector<Atom>::iterator at;
    std::vector<Coordinate>::iterator ct=Molecules()->back().Frames().back().XYZ().begin();
    for ( at=Molecules()->back().Atoms().begin(),i=1;at!=Molecules()->back().Atoms().end();++at,++ct,i++ )
    {
      std::getline ( *m_file,line );
      StringTokenizer token ( line );
      if ( token.size() < 4 )
      {
        throw kryomol::Exception ( "Error parsering mdl V2000 file" );
      }
      ct->x() =std::atof ( token.at(0).c_str() );
      ct->y() =std::atof ( token.at(1).c_str() );
      ct->z() =std::atof ( token.at(2).c_str() );
      std::stringstream st;
      st << i;
      at->SetSymbol ( token.at(3) );

    }


  }
  return true;

}

