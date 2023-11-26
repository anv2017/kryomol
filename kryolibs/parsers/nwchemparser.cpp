/*****************************************************************************************
                            nwchemparser.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <clocale>
#include "exception.h"
#include "molecule.h"
#include "stringtools.h"

#include "nwchemparser.h"

using namespace kryomol;

NwChemParser::NwChemParser ( const char* file )
    : Parser ( file )
{
  m_nsteps=0;
}

NwChemParser::NwChemParser ( std::istream* stream ) : Parser ( stream )
{}



NwChemParser::~NwChemParser()
{}

bool NwChemParser::ParseFile ( std::streampos pos /*=0*/ )
{
  m_file->clear();
  m_file->seekg ( pos,std::ios::beg );
  GetOptimizedGeometry();
  if ( Molecules()->empty() )
  {
    m_file->clear();
    m_file->seekg ( pos,std::ios::beg );
    GetSinglePointGeometry();
  }
  return true;

}

bool NwChemParser::GetSinglePointGeometry()
{
  std::string line;
  Molecules()->resize ( 1 );
  Molecules()->back().Frames().push_back ( Frame ( &Molecules()->back() ) );
  while ( std::getline ( *m_file,line ) )
  {
    if ( line.find ( "Geometry" ) != std::string::npos )
    {

      for ( int i=0;i<6;++i )
        std::getline ( *m_file,line );
      while ( std::getline ( *m_file,line ) )
      {
        StringTokenizer token ( line," \t\r" );
        if ( token.size() < 6 )
        {
          return true;
        }

        Atom atom ( token.at ( 1 ) );;
        Coordinate c;
        c.x() =kryomol::atof ( token.at ( 3 ) );
        c.y() =kryomol::atof ( token.at ( 4 ) );;
        c.z() =kryomol::atof ( token.at ( 5 ) );
        Molecules()->back().Atoms().push_back ( atom );
        Molecules()->back().Frames().back().XYZ().push_back ( c );
      }
    }
  }
  return true;
}

bool NwChemParser::GetOptimizedGeometry()
{
  std::streampos orgpos=m_file->tellg();
  std::string line;
  std::streampos pos=0;
  //Find the last step
  while ( std::getline ( *m_file,line ) )
  {
    if ( line.find ( "Step" ) != std::string::npos )
      pos=m_file->tellg();

    if ( line.find ( "Optimization converged" ) != std::string::npos ) return true;
  }
  //work this on windows
  if ( pos > 0 )
  {
    m_file->clear();
    m_file->seekg ( pos );


    if ( !GetGeometryBlock() ) return false;
    // GetEnergies();
  }

  if ( Molecules()->empty() )
  {
    m_file->clear();
    m_file->seekg ( orgpos );
    pos=0;
    while ( std::getline ( *m_file,line ) )
    {
      if ( line.find ( "Atom information" ) != std::string::npos )
        pos=m_file->tellg();
    }
    if ( pos > 0 )
   {
    m_file->clear();
    m_file->seekg ( pos );

    while ( std::getline ( *m_file,line ) )
    {
      if ( line.find ( "Atom information" ) != std::string::npos )
      {
        std::getline ( *m_file,line );
        std::getline ( *m_file,line );
        Molecules()->resize ( 1 );
        Molecules()->back().Frames().push_back ( Frame ( &Molecules()->back() ) );
        while ( std::getline ( *m_file,line ) )
        {
          if ( line.find ( "---" ) !=std::string::npos ) break;
          StringTokenizer token ( line," \t\r" );
          Atom atom ( token.at ( 0 ) );
          Coordinate c;
          c.x() =kryomol::atof ( token.at ( 2 ) );
          c.y() =kryomol::atof ( token.at ( 3 ) );
          c.z() =kryomol::atof ( token.at ( 4 ) );
          c/=1.889725989;

          Molecules()->back().Atoms().push_back ( atom );
          Molecules()->back().Frames().back().XYZ().push_back ( c );
        }
        return true;
      }
    }
  }
  else return false;
  }
  return true;
}

bool NwChemParser::GetGeometryBlock()
{
  //skip three lines
  std::string line;
  for ( int i=0;i<3;i++ )
  {
    std::getline ( *m_file,line );
  }
  Molecules()->resize ( 1 );
  Molecules()->back().Frames().push_back ( Frame ( &Molecules()->back() ) );

  while ( std::getline ( *m_file,line ) )
  {
    StringTokenizer token ( line," \t\r" );
    if ( token.size() )
    {
      Atom atom ( token.at ( 1 ) );
      Coordinate c;
      c.x() =std::atof ( token.at ( 2 ).c_str() );
      c.y() =std::atof ( token.at ( 3 ).c_str() );
      c.z() =std::atof ( token.at ( 4 ).c_str() );
      c/=1.889725989;

      Molecules()->back().Atoms().push_back ( atom );
      Molecules()->back().Frames().back().XYZ().push_back ( c );
    }
    else break;

  }

  return true;
}

void NwChemParser::ParseCouplingConstants ( std::vector<QuantumCoupling>& c )
{
  const char* locale=std::setlocale(LC_NUMERIC,"C");
  std::setlocale(LC_NUMERIC,locale);
  m_file->clear();
  m_file->seekg ( 0,std::ios::beg );
  std::string line;
  while ( std::getline ( *m_file,line ) )
  {
    if ( line.find ( "Indirect Spin-Spin Tensors" ) != std::string::npos )
    {
      while ( std::getline ( *m_file,line ) )
      {
        size_t i,j;
        if ( line.find ( "Atom" ) != std::string::npos )
        {

          StringTokenizer tok ( line," :\t\r" );
          i = atoi ( tok.at ( 1 ).c_str() )-1;
          j = atoi ( tok.at ( 5 ).c_str() )-1;
          c.push_back ( QuantumCoupling ( i ,j ) );

        }

        if ( line.find ( "Fermi Contact" ) )
        {
          for ( int i=0;i<4;++i )
            std::getline ( *m_file,line );

          if ( line.find ( "Isotropic" ) != std::string::npos )
          {
            StringTokenizer tok ( line," \t\r" );
            c.back().SetFermiContact ( kryomol::atof ( tok.back() ) );

          }
          else throw kryomol::Exception ( "Error parsering fermi contact contributions" );

        }

        if ( line.find ( "Spin-Dipole" ) )
        {
          for ( int i=0;i<4;++i )
            std::getline ( *m_file,line );

          if ( line.find ( "Isotropic" ) != std::string::npos )
          {
            StringTokenizer tok ( line," \t\r" );
            c.back().SetSpinDipole ( kryomol::atof ( tok.back() ) );

          }
          else throw kryomol::Exception ( "Error parsering spin-dipole contributions" );

        }

        if ( line.find ( "Cross Term" ) )
        {
          for ( int i=0;i<4;++i )
            std::getline ( *m_file,line );

          if ( line.find ( "Isotropic" ) != std::string::npos )
          {
            StringTokenizer tok ( line," \t\r" );
            c.back().SetFCSDCrossTerm ( kryomol::atof ( tok.back() ) );

          }
          else throw kryomol::Exception ( "Error parsering Fermi-Contact Spin Dipole cross term contributions" );

        }

        if ( line.find ( "Paramagnetic" ) )
        {
          for ( int i=0;i<4;++i )
            std::getline ( *m_file,line );

          if ( line.find ( "Isotropic" ) != std::string::npos )
          {
            StringTokenizer tok ( line," \t\r" );
            c.back().SetParamagneticSpinOrbit ( kryomol::atof ( tok.back() ) );

          }
          else throw kryomol::Exception ( "Error parsering diamagnetic Spin-Orbit contributions" );

        }


        if ( line.find ( "Diamagnetic" ) )
        {
          for ( int i=0;i<4;++i )
            std::getline ( *m_file,line );

          if ( line.find ( "Isotropic" ) != std::string::npos )
          {
            StringTokenizer tok ( line," \t\r" );
            c.back().SetDiamagneticSpinOrbit ( kryomol::atof ( tok.back() ) );

          }
          else throw kryomol::Exception ( "Error parsering diamagnetic Spin-Orbit contributions" );

        }

        if ( line.find ( "Total Spin-Spin" ) != std::string::npos && ( i < j ) )
        {
          for ( int i=0;i<4;++i )
            std::getline ( *m_file,line );

          if ( line.find ( "Isotropic" ) != std::string::npos )
          {
            StringTokenizer tok ( line," \t\r" );
            c.back().SetValue( kryomol::atof ( tok.back() ) );

          }
          else throw kryomol::Exception ( "Error parsering total spin-spin couplings contributions" );

        }
      }
    }

  }

  std::setlocale(LC_NUMERIC,locale);
}

// std::vector<Quantum::jobheader>& NwChemParser::GetJobs()
// {
//   //lets get the lines beginning with a #
//   std::string line;
//   std::streampos pos;
//   Quantum::jobtype tjob;
//   while ( std::getline ( *m_file,line ) )
//   {
//     bool bfound=false;
//
//     if ( line.find ( "NWChem Geometry Optimization",0 ) != std::string::npos )
//     {
//       tjob=Quantum::opt;
//       bfound=true;
//     }
//
//     if ( line.find ( " NWChem Nuclear Hessian and Frequency Analysis",0 ) != std::string::npos )
//     {
//       tjob=Quantum::freq;
//       bfound=true;
//     }
//
//     if ( bfound )
//     {
//       pos=m_file->tellg();
//       m_jobpos.push_back ( Quantum::jobheader ( tjob,pos ) );
//     }
//
//   }
//   if  ( m_jobpos.empty() )
//   {
//     m_jobpos.push_back( Quantum::jobheader( Quantum::dens,0 ));
//   }
//
//   return m_jobpos;
// }
