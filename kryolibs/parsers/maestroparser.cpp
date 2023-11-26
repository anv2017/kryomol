/*****************************************************************************************
                            maestroparser.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/


#include "maestroparser.h"
#include "molecule.h"
#include "stringtools.h"
#include "exception.h"

const double rchsp3=1.100;
const double rchsp2=1.090;
const double rchsp=1.080;
const double tetsin=0.8164966;
const double tetcos=0.5773503;

using namespace kryomol;

MaestroParser::MaestroParser ( const char* file )
    : Parser ( file )
{}

MaestroParser::MaestroParser ( std::istream* stream ) : Parser ( stream )
{
}

MaestroParser::~MaestroParser()
{}

bool MaestroParser::ParseFile ( std::streampos pos )
{
  m_file->clear();
  m_file->seekg ( pos );
  GetGeometries();

  return true;

}

void MaestroParser::GetGeometries()
{
  std::string line;
  Molecules()->push_back ( Molecule() );
  while ( std::getline ( *m_file,line ) )
  {
    if ( line.find ( "f_m_ct" ) != std::string::npos )
    {
      GetBlock ( true );
    }
    if ( line.find ( "p_m_ct" ) != std::string::npos )
    {
      GetBlock ( false );
    }
  }

  if ( Molecules()->back().Frames().back().XYZ().size() > Molecules()->back().Atoms().size() )
  {
    int newh=Molecules()->back().Frames().back().XYZ().size() - Molecules()->back().Atoms().size();
    for ( int i=0;i<newh;++i )
    {
      Molecules()->back().Atoms().push_back ( Atom ( 1 ) );
    }
  }
//  Molecules()->back().MoveToCentroid();

}


void MaestroParser::GetBlock ( bool full )
{
  int counter=0;
  bool hasenergy=false;
  std::string line;

  Molecules()->back().Frames().push_back ( Frame ( &Molecules()->back() ) );
  Frame& frame=Molecules()->back().Frames().back();
  while ( std::getline ( *m_file,line ) )
  {
    if ( !hasenergy )
      counter++;
    if ( line.find ( "r_mmod_Potential_Energy" ) != std::string::npos ) hasenergy=true;
    if ( line.find ( ":::" ) != std::string::npos ) break;
  }
  if ( hasenergy )
  {
    for ( int i=0;i<counter;i++ )
    {
      std::getline ( *m_file,line );
    }

    StringTokenizer token ( line," \t\r" );
    frame.PotentialEnergy() = ( std::atof ( token.back().c_str() ) *0.2390 );

  }
  //find the atoms
  while ( std::getline ( *m_file,line ) )
  {
    if ( line.find ( "m_atom" ) != std::string::npos )
    {
      StringTokenizer token ( line," \t\r[]{}" );

      if ( Molecules()->back().Atoms().empty() )
        Molecules()->back().Atoms().resize ( atoi ( token.back().c_str() ) );

      frame.XYZ().resize ( atoi ( token.back().c_str() ) );
      goto lb;
    }

  }
lb:
  while ( std::getline ( *m_file,line ) )
  {
    if ( line.find ( ":::" ) != std::string::npos ) break;

  }

  //OK lets do the parsering
  while ( std::getline ( *m_file,line ) )
  {
    if ( line.find ( ":::" ) != std::string::npos ) break;
    StringTokenizer token ( line," \t\r" );
    Atom& atom=Molecules()->back().Atoms().at ( atoi ( token.at ( 0 ).c_str() )-1 );
    Coordinate& c=frame.XYZ().at ( atoi ( token.at ( 0 ).c_str() )-1 );
    if ( !full )
    {
      c.x() =std::atof ( token.at ( 1 ).c_str() );
      c.y() =std::atof ( token.at ( 2 ).c_str() );
      c.z() =std::atof ( token.at ( 3 ).c_str() );
    }
    else
    {
      int z=AtomType ( atoi ( token.at ( 1 ).c_str() ) );
      if ( z == 0 ) throw Exception("Maestro atom type not recognized");
      if ( z == -1 )
      {
        m_uatoms.push_back ( std::pair<size_t,int> ( atoi ( token.at ( 0 ).c_str() )-1,atoi ( token.at ( 1 ).c_str() ) ) );
        z=6;
      }
      atom.SetZ ( z );
      c.x() =std::atof ( token.at ( 2 ).c_str() );
      c.y() =std::atof ( token.at ( 3 ).c_str() );
      c.z() =std::atof ( token.at ( 4 ).c_str() );
    }
    if ( line.find ( ":::" ) != std::string::npos ) break;
  }

  if ( !m_uatoms.empty() )
  {
    Molecules()->back().SetBonds();
    TransformUnitedAtoms();

  }

}

int MaestroParser::AtomType ( int type )
{
  switch ( type )
  {
    case 1: //carbon sp
    case 2:  //carbon sp2
    case 3:  //carbon sp3
      return 6;
      break;
    case 4:  //UA CH-sp3
    case 5:  //UA CH2-sp3
    case 6:  //UA CH3-sp3
    case 7:  //UA CH-sp2
    case 8:  //UA CH2-sp2
    case 9:  //UA CH-sp
      //throw error("Not support yet for united atoms");
      return -1;
    case 10: //C-
    case 11: //C+
    case 12: //C
    case 14: //any carbon
      return 6;  //carbon sp2
    case 15: //O single bond
    case 16: //O single bond
      return 8;
    case 17:
      // throw error("Not support yet for united atoms");
      return -1;
    case 18: //O-
      return 8;
    case 19: //H2O united ignore
      return 0;
    case 20: //R2=O+
    case 21: //R3O+
    case 23: //any oxygen
      return 8;
    case 24: //N sp
    case 25: //N sp2
    case 26: //N sp3
      return 7;
    case 27:  //UA NH-sp3
    case 28:  //UA NH2-sp3
    case 29:  //UA NH-sp2
    case 30:  //UA NH2-sp2
      // throw error("Not support yet for united atoms");
      return -1;
    case 31:  //N+-sp2
    case 32:  //N+-sp3
      return 7;
    case 33:  //UA NH+-sp3
    case 34:  //UA NH2+-sp3
    case 35:  //UA NH3+-sp3
    case 36: //UA NH+-sp2
    case 37: //UA NH2+-sp2
      // throw error("Not support yet for united atoms");
      return -1;
    case 38: //N- sp3
    case 39: //N- sp2
    case 40: //Any nitrogen
      return 7;
    case 41: //H-Electroneut (C,S)
    case 42: //H-O
    case 43: //H-N
    case 44: //H-cation
    case 45: //H-Anion
    case 48: //Any hydrogen
      return 1;
    case 49: //sulfur
      return 16;
    case 50: //United atom SH
//   throw error("Not support yet for united atoms");
      return -1;
    case 51: //S-
    case 52: //any sulphur
      return 16;
    case 53: //P
      return 15;
    case 54: //B sp2
    case 55: //B sp3
      return 5;
    case 56: //F
      return 9;
    case 57: //Cl
      return 17;
    case 58: //Br
      return 35;
    case 59: //I
      return 53;
    case 60: //Si
      return 14;
    case 61: //dummy atom for FEP
    case 62: //Special atom to be defined
    case 63: //Lp
    case 64: //Any atom
      return 0;
    case 65: //Li+
      return 3;
    case 66: //Na++
      return 11;
    case 67: //K+
      return 19;
    case 68://Rb+
      return 37;
    case 69: //Cs+
      return 55;
    case 70: //Ca++
      return 20;
    case 71: //Ba++
      return 56;
    case 72: //Mg++
      return 12;
    case 73: //Mn+2
    case 74: //Mn+3
    case 75: //Mn+4
    case 76: //Mn+5
    case 77: //Mn+6
    case 78: //Mn+7
      return 25;
    case 79: //Fe+2
    case 80: //Fe+3
      return 26;
    case 81: //Co+2
    case 82: //Co+3
      return 27;
    case 83: //Ni+2
    case 84: //Ni+3
      return 28;
    case 85: //Cu+
    case 86: //Cu+2
      return 29;
    case 87: //Zn+2
      return 30;
    case 88: //Mo+3
    case 89: //Mo+4
    case 90: //Mo+5
    case 91: //Mo+6
      return 42;
    case 100: //S+
    case 101: //S sp2
      return 16;
    case 102: //Cl-
      return 17;
    case 109: // tetrahedral
    case 110: //S octahedral
        return 16;
    case 111:  //posphorous cation
        return 15;
    case 112: //selenium
         return 34;
    case 113: //sulphur hexavalent tetrahedral
    case 114: //sulphur sulfide anion
         return 16;
    default:
      return 0;
  }
  return 0;
}

void MaestroParser::TransformUnitedAtoms()
{
  std::vector< std::pair<size_t,int> >::iterator it;
  for ( it=m_uatoms.begin();it!=m_uatoms.end();++it )
  {
    
    TransformAtom ( *it );
    std::cout << "Hydrogens have been added to united atom " << it->first+1 << std::endl;
  }
  std::vector<Coordinate>::const_iterator ct;
  for ( ct=m_newh.begin();ct!=m_newh.end();++ct )
  {
    
    Molecules()->back().Frames().back().XYZ().push_back ( *ct );
  }
  m_newh.clear();
}

void MaestroParser::TransformAtom ( const std::pair<size_t,int>& atom )
{
  switch ( atom.second )
  {
    case 4:
      AddHToCR3 ( atom.first );
      break;
    case 5:
      AddHToCR2 ( atom.first );
      break;
    case 6:
      AddHToCR ( atom.first );
      break;
    case 7:
      AddHToCsp2 ( atom.first );
      break;
    case 8:
      AddH2ToCsp2 ( atom.first );
      break;
    case 9:
      AddHToCsp ( atom.first );
      break;
    default:
      throw kryomol::Exception ( "United atom not supported" );
      break;
  }
}

void MaestroParser::AddHToCR3 ( size_t atom )
{
  const Coordinate& c=Molecules()->back().Frames().back().XYZ().at ( atom );
  std::vector<size_t> nh=Molecules()->back().Neighbours ( atom );
  std::vector<size_t>::const_iterator it;
  const Coordinate& c1=Molecules()->back().Frames().back().XYZ().at ( nh.at ( 0 ) );
  const Coordinate& c2=Molecules()->back().Frames().back().XYZ().at ( nh.at ( 1 ) );
  const Coordinate& c3=Molecules()->back().Frames().back().XYZ().at ( nh.at ( 2 ) );
  Coordinate v1=c1-c;
  Coordinate v2=c2-c;
  Coordinate v3=c3-c;
  v1/=v1.Norm();
  v2/=v2.Norm();
  v3/=v3.Norm();
  Coordinate s=v1+v2+v3;
  s*=-1;
  s/=s.Norm();
  s*=rchsp3;
  m_newh.push_back ( s+c );
}

void MaestroParser::AddHToCR2 ( size_t atom )
{
  const Coordinate& c=Molecules()->back().Frames().back().XYZ().at ( atom );
  std::vector<size_t> nh=Molecules()->back().Neighbours ( atom );
  std::vector<size_t>::const_iterator it;
  const Coordinate& c1=Molecules()->back().Frames().back().XYZ().at ( nh.at ( 0 ) );
  const Coordinate& c2=Molecules()->back().Frames().back().XYZ().at ( nh.at ( 1 ) );
  Coordinate v1=c1-c;
  Coordinate v2=c2-c;
  v1/=v1.Norm();
  v2/=v2.Norm();
  v1=v1+v2;
  v2=v1^v2;
  v1=v1/v1.Norm();
  v1*=-1;
  v2=v2/v2.Norm();
  Coordinate s=c+v1*tetcos+v2*tetsin;
  Coordinate s1=c+v1*tetcos-v2*tetsin;
  s*=rchsp3;
  s1*rchsp3;
  m_newh.push_back ( s );
  m_newh.push_back ( s1 );

}

//Add methyl group
void MaestroParser::AddHToCR ( size_t atom )
{
  const Coordinate& c=Molecules()->back().Frames().back().XYZ().at ( atom );
  std::vector<size_t> nh=Molecules()->back().Neighbours ( atom );
  std::vector<size_t>::const_iterator it;
  const Coordinate& c1=Molecules()->back().Frames().back().XYZ().at ( nh.at ( 0 ) );

  double prjz=-rchsp3*sin ( ( 180-109.47122 ) /M_PI );
  double prjx=-rchsp3*cos ( ( 180-109.47122 ) /M_PI );
  Coordinate h1=Coordinate ( prjx,0,prjz ) +c;
  Coordinate h2=Coordinate ( -prjx*0.5,0.8660*prjx,prjz ) +c;
  Coordinate h3=Coordinate ( -prjx*0.5,-0.8660*prjx,prjz ) +c;

  Coordinate v1=c1-c;
  Coordinate z ( 0,0,1 );
  Coordinate v2=v1^z;
  double angle=Coordinate::Angle ( v1,z );
#ifdef __GNUC__
#warning check that this is still correct
#endif
  //c,v2+c
  h1=Coordinate::RotAroundAxis ( c,v2+c,h1,M_PI-angle );
  h2=Coordinate::RotAroundAxis ( c,v2+c,h2,M_PI-angle );
  h3=Coordinate::RotAroundAxis ( c,v2+c,h3,M_PI-angle );
  m_newh.push_back ( h1 );
  m_newh.push_back ( h2 );
  m_newh.push_back ( h3 );


}

void MaestroParser::AddHToCsp ( size_t atom )
{
  const Coordinate& c=Molecules()->back().Frames().back().XYZ().at ( atom );
  std::vector<size_t> nh=Molecules()->back().Neighbours ( atom );
  std::vector<size_t>::const_iterator it;
  const Coordinate& c1=Molecules()->back().Frames().back().XYZ().at ( nh.at ( 0 ) );

  Coordinate v1=c1-c;
  v1*=-1;
  v1/=v1.Norm();
  v1*=rchsp;
  Coordinate s1= c+v1 ;
  m_newh.push_back ( s1 );

}

void MaestroParser::AddH2ToCsp2 ( size_t atom )
{
  const Coordinate& c=Molecules()->back().Frames().back().XYZ().at ( atom );
  std::vector<size_t> nh=Molecules()->back().Neighbours ( atom );
  std::vector<size_t>::const_iterator it;
  const Coordinate& c1=Molecules()->back().Frames().back().XYZ().at ( nh.at ( 0 ) );

  std::vector<size_t> nh1=Molecules()->back().Neighbours ( nh.at ( 0 ) );
  size_t tatom = nh1.at ( 0 );
  if ( tatom == atom ) tatom=nh1.at ( 1 );
  Coordinate c2=Molecules()->back().Frames().back().XYZ().at ( tatom );
  Coordinate v1=c2-c1;
  Coordinate v2=c-c1;
  Coordinate v3=v1^v2;
  Coordinate v=v3^v2;
  v/=v.Norm();
  v2/=v2.Norm();
  Coordinate s1=c+v*rchsp2*cos ( 30*M_PI/180. ) +v2*rchsp2*0.5;
  Coordinate s2=c-v*rchsp2*cos ( 30*M_PI/180. ) +v2*rchsp2*0.5;
  m_newh.push_back ( s1 );
  m_newh.push_back ( s2 );
}

void MaestroParser::AddHToCsp2 ( size_t atom )
{
  const Coordinate& c=Molecules()->back().Frames().back().XYZ().at ( atom );
  std::vector<size_t> nh=Molecules()->back().Neighbours ( atom );
  std::vector<size_t>::const_iterator it;
  const Coordinate& c1=Molecules()->back().Frames().back().XYZ().at ( nh.at ( 0 ) );
  const Coordinate& c2=Molecules()->back().Frames().back().XYZ().at ( nh.at ( 1 ) );

  Coordinate v1=c2-c;
  Coordinate v2=c1-c;
  Coordinate v=v1+v2;
  v/=v.Norm();

  Coordinate s1=c-v*rchsp2;
  m_newh.push_back ( s1 );

}

