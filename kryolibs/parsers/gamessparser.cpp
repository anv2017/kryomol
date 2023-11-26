/*****************************************************************************************
                            gamesparser.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <cstdlib>
#include <clocale>

#include "exception.h"
#include "stringtools.h"
#include "gamessparser.h"
#include "molecule.h"


const double bohrtoangs=0.529189379;

using namespace kryomol;

GamessParser::GamessParser ( const char* file ) : Parser ( file )
{}

GamessParser::GamessParser ( std::istream* stream ) : Parser ( stream )
{}


GamessParser::~GamessParser()
{}

bool GamessParser::ParseFile ( std::streampos pos )
{
  m_file->clear();
  m_file->seekg ( pos,std::ios::beg );
  GetGeometry();

  return true;

}

void GamessParser::GetGeometry()
{
  std::string line;
  std::streampos pos =0;
  bool bfound=false;
  bool bbohr=true;

  while ( std::getline ( *m_file,line ) )
  {

    if ( line.find ( "COORDINATES",0 ) != std::string::npos && line.find ( "INTERNAL" ) == std::string::npos ) //this is for optimizations
    {
      short skiplines=2;
      if ( line.find ( "BOHR",0 ) != std::string::npos )
      {
        bbohr=true;
        skiplines=1;
      }
      bfound=true;
      for ( int i=0;i<skiplines;i++ )
      {
        std::getline ( *m_file,line );

      }
      pos= m_file->tellg() ;

    }
  }

  if ( !bfound ) throw kryomol::Exception ( "Geometry not found" );

  m_file->clear();
  m_file->seekg ( pos );
  Molecules()->resize ( 1 );
  Molecules()->back().Frames().push_back ( Frame ( &Molecules()->back() ) );

  while ( std::getline ( *m_file,line ) )
  {
    StringTokenizer token ( line," \t\r" );
    size_t size=token.size();

    if ( size >= 5 )
    {

      Atom atom ( atoi ( token.at ( 1 ).c_str() ) ) ;
      Coordinate c;
      c.x() =std::atof ( token.at ( size-3 ).c_str() );
      c.y() =std::atof ( token.at ( size-2 ).c_str() );
      c.z() =std::atof ( token.at ( size-1 ).c_str() );

      if ( bbohr ) c*=bohrtoangs;

      Molecules()->back().Atoms().push_back ( atom );
      Molecules()->back().Frames().back().XYZ().push_back ( c );



    }

    else
      break;
  }


}


void GamessParser::ParseChemicalShifts()
{
  const char* locale=std::setlocale(LC_NUMERIC,"C");

  m_file->clear();
  m_file->seekg ( 0 );
  std::string line;
  while ( std::getline ( *m_file,line ) )
  {
    if ( line.find ( "CHEMICAL SHIELDING TENSOR" ) != std::string::npos )
    {

      while ( std::getline ( *m_file,line ) )
      {

        if ( line.find ( "(" ) != std::string::npos )
        {
          std::getline ( *m_file,line );
          if ( !ExtractTensorBlock() ) return;
        }
      }

    }
  }
  std::setlocale(LC_NUMERIC,locale);
}

bool  GamessParser::ExtractTensorBlock()
{
  Frame& frame=Molecules()->back().Frames().back();
  D2Array<double> tensor ( 3,3 ) ;
  size_t row=0;
  std::string line;
  for ( int i=0;i<3;++i )
  {
    StringTokenizer  tok ( line," \t\r" );
    
    if ( tok.empty() || !kryomol::isinteger( tok.at(0) ) ) return false;
    StringTokenizer::reverse_iterator it;
    int counter=0;

    for ( it=tok.rbegin();it!=tok.rend();++it,++counter )
    {
      int col=2-counter;
      tensor( row,col ) = kryomol::atof(*it);
    }
    ++row;
  }
  
  frame.CShiftTensors().push_back(tensor);
  return true;
}

// std::vector<Quantum::jobheader>& GamessParser::GetJobs()
// {
//   std::string line;
//   Quantum::jobtype tjob=Quantum::dens;
//   while ( std::getline ( *m_file,line ) )
//   {
//     if ( line.find ( "$CONTRL OPTIONS" ) != std::string::npos )
//       break;
//   }
//
//   std::getline ( *m_file,line );
//   std::getline ( *m_file,line );
//
//   if ( line.find ( "RUNTYP=HESSIAN" ) != std::string::npos )
//     tjob=Quantum::freq;
//   if ( line.find ( "RUNTYP=OPTIMIZE" ) != std::string::npos )
//     tjob=Quantum::opt;
//
//   m_jobpos.push_back ( Quantum::jobheader ( tjob,0 ) );
//
//   return m_jobpos;
//
// }

