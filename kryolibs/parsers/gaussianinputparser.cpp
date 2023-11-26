/*****************************************************************************************
                            gaussianinputparser.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <sstream>

#include "molecule.h"
#include "stringtools.h"
#include "gaussianinputparser.h"

using namespace kryomol;

GaussianInputParser::GaussianInputParser ( const char* file ) : Parser ( file )
{}

GaussianInputParser::GaussianInputParser ( std::istream* stream ) :Parser ( stream )
{}

GaussianInputParser::~GaussianInputParser()
{}


bool GaussianInputParser::ParseFile ( std::streampos pos )
{

  m_file->clear();
  m_file->seekg ( pos,std::ios::beg );

  GetGeometry();

  return true;
}


void GaussianInputParser::GetGeometry()
{
  std::string line;
  Molecules()->push_back ( Molecule() );
  Molecules()->back().Frames().push_back ( Frame ( &Molecules()->back() ) );
  while ( std::getline ( *m_file,line ) )
  {
    StringTokenizer token ( line," \t\r," );
    if ( token.size() == 2 )
    {

      if ( isinteger ( token.at ( 0 ) ) && isinteger ( token.at ( 1 ) ) )
      {
        while ( std::getline ( *m_file,line ) )
        {

          StringTokenizer satom ( line," \t\r," );
          if ( satom.size() < 4 ) break;
          if ( satom.at ( 0 ) !="Lp" ) //skip lone pairs
          {
            Atom atom ( satom.front() );
            Coordinate c;
            StringTokenizer::reverse_iterator it2=satom.rbegin();
            c.z() = std::atof ( ( it2++ )->c_str() );
            c.y() = std::atof ( ( it2++ )->c_str() );
            c.x() = std::atof ( it2->c_str() );
            Molecules()->back().Atoms().push_back ( atom );
            Molecules()->back().Frames().back().XYZ().push_back ( c );
          }


        }


      }
    }
  }

}
