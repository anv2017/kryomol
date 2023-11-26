/*****************************************************************************************
                            acesparser.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <sstream>

#include "acesparser.h"
#include "stringtools.h"
#include "molecule.h"


const double bohrtoangs=0.529189379;

using namespace kryomol;

AcesParser::AcesParser(const char* file) : Parser(file)
{}

AcesParser::AcesParser(std::istream* stream) : Parser(stream)
{}

AcesParser::~AcesParser()
{}

bool AcesParser::ParseFile(std::streampos pos)
{
  m_file->clear();
  m_file->seekg(0/*pos*/,std::ios::beg);
  m_level=GetLevel();
  std::cerr << "passed get level" << std::endl;
  bool b = GetGeometry();

  return b;
}

bool AcesParser::GetGeometry()
{

  std::string line;
  bool bfound=false;
  std::vector<std::streampos> pos;
  while(std::getline(*m_file,line))
  {
    StringTokenizer tok(line," ");
    if( tok.size() >= 2)
    {
      if( tok[0]=="Z-matrix" && tok[1]=="Atomic")
      {
        std::cerr << "found" << std::endl;
        bfound=true;
        pos.push_back(m_file->tellg());

      }
    }
  }

  std::cerr << "In Aces Parser" << std::endl;

  if(!bfound) return false;

  if ( !pos.empty() )
  {
      std::cout << "pushing molecule " << std::endl;
          Molecules()->push_back ( Molecule() );
  }


  for(std::vector<std::streampos>::iterator pt=pos.begin();pt!=pos.end();++pt)
  {
    m_file->clear();
    m_file->seekg(*pt,std::ios::beg);
    std::getline(*m_file,line);
    std::getline(*m_file,line);

    Molecules()->back().Frames().push_back ( Frame ( &Molecules()->back() ) );

    while(std::getline(*m_file,line))
    {
      StringTokenizer token(line," ");
      size_t size=token.size();
      if( size == 5  )
      {
        if ( token.at(0) != "X" )
        {
        std::string symbol = token.at(0);
        //Aces writes all atoms as uppercase
        if(symbol.size()> 1)
        {
          symbol.at(1)=std::tolower(symbol.at(1),std::cout.getloc());
        }
        if ( pt == pos.begin() )
        {
          Molecules()->back().Atoms().push_back( Atom (symbol) );
        }
        Coordinate c;
        c.x() = std::atof(token[size-3].c_str())*bohrtoangs;
        c.y() = std::atof(token[size-2].c_str())*bohrtoangs;
        c.z() = std::atof(token[size-1].c_str())*bohrtoangs;
        Molecules()->back().Frames().back().XYZ().push_back ( c );
        }
      }

      else
        break;
    }

    while(std::getline(*m_file,line))
    {
      bool bfound=false;
      StringTokenizer token(line," ");
      switch(m_level)
      {
      case HF:
        if(line.find("E(SCF)")!=std::string::npos)
        {
          StringTokenizer token(line," ");
          Molecules()->back().Frames().back().SetEnergy ( std::atof ( token.at ( 1 ).c_str() ), "SCF" );//.c_str() ),"SCF" );
        }
        break;
      case CCSD:
         if ( line.find("Total CCSD energy") != std::string::npos )
         {
                StringTokenizer token(line, " \t\r");
                Molecules()->back().Frames().back().SetEnergy ( std::atof ( token.back().c_str()) );
         }
         break;
      case CCSDT:
        if( token.size() == 3 )
        {
          if( token.at(0)=="CCSD(T)" &&  ( token.at(1)=="=" || token.at(1)=="energy" ) )
            {
            std::cout << "energy(ccsdt)=" << std::atof(token.back().c_str() );
            Molecules()->back().Frames().back().SetEnergy ( std::atof ( token.back().c_str()) );
        }
        }
        break;
      default:
        break;

      }
      if(line.find("Minimum force",0)!=std::string::npos)
      {
        StringTokenizer token(line," ");
        Molecules()->back().Frames().back().SetRMSForce ( std::atof ( token.back().c_str() ) );
        if(token.size() > 2 )
            Molecules()->back().Frames().back().SetMaximumForce ( std::atof ( token.at ( 2 ).c_str() ) );
        bfound=true;
      }
      if(bfound) break;
    }


  }
  Molecules()->back().SetBonds();
  return true;
}


QuantumLevel AcesParser::GetLevel()
{
  std::string line;
  while(std::getline(*m_file,line))
  {
    if(line.find("CALCLEVEL",0)!=std::string::npos)
      break;
  }

  StringTokenizer token(line," []");

  StringTokenizer::reverse_iterator it=token.rbegin();
  it++;

  std::string level;
  if(it!=token.rend())
    level=(*it);


  //leave in a good state
  m_file->clear();
  m_file->seekg(0,std::ios::beg);

  std::cerr << "leaving get level=" << level << std::endl;
  if(level=="22") return CCSDT;
  if(level=="10") return CCSD;
  if(level=="0") return HF;
  if(level=="1") return MP2;

  return HF;

}


bool AcesParser::ParseFrequencies(std::streampos pos)
{
  m_file->clear();
  m_file->seekg(pos,std::ios::beg);

  Molecule& molecule=Molecules()->back();
  molecule.Frames().back().AllocateVectors();

  m_counter=0;

  GetFrequencies();

  GetVectors();

  return true;
}


void AcesParser::GetFrequencies()
{
  Molecules()->back().Frames().back().AllocateVectors();
  std::vector<Frequency>& frequencies = Molecules()->back().Frames().back().GetFrequencies();

  std::string line;
  while(std::getline(*m_file,line))
  {
    if(line.find("Representation       Frequency ") != std::string::npos)
      break;
  }
  for(int i=0;i<3;i++)
    std::getline(*m_file,line);

  while(std::getline(*m_file,line))
  {
    StringTokenizer token(line," \t");
    if ( token.size() < 4 ) break;
    if( token.back() == "VIBRATION")
    {
        frequencies.push_back( Frequency( std::atof(token.at(1).c_str()), std::atof(token.at(2).c_str()),0));
    }
  }
}


void AcesParser::GetVectors()
{
  std::string line;
  while(std::getline(*m_file,line))
  {
   // if(line.find("X       Y       Z") != std::string::npos)
    if ( line.find("VIBRATION") != std::string::npos )
    {

      StringTokenizer token(line," \t\r");
      if( !GetVectorBlock(token.size()) ) break;
    }

  }

}

bool AcesParser::GetVectorBlock(int nmodes)
{
  int atom=0;
  std::string line;

  std::getline(*m_file,line);
  if ( line.find("X       Y       Z") != std::string::npos ) std::getline(*m_file,line);
do
{
    if (line.find("Gradient vector in normal") != std::string::npos) return false;
    //There is a bug in ACES formatting so insert a space after 12th character
    if(line.size() > 13 ) line.insert(13," ");
    StringTokenizer token(line," \t");
    if( token.size() < 4) break;
    std::vector< std::vector<Coordinate> >& coordinates = Molecules()->back().GetMode();
    for(int j=0;j<nmodes;j++)
    {
      coordinates.at(j+m_counter).at(atom).x() = std::atof ( token[1+3*j].c_str() );
      coordinates.at(j+m_counter).at(atom).y() = std::atof ( token[2+3*j].c_str() );
      coordinates.at(j+m_counter).at(atom).z() = std::atof ( token[3+3*j].c_str() );
    }
    atom++;

  } while (std::getline(*m_file,line) );

  m_counter+=nmodes;
  return true;
}

std::vector<JobHeader>& AcesParser::Jobs()
{
  //lets get the lines beginning with a #
  m_file->clear();
  m_file->seekg(0,std::ios::beg);
  std::string line;
  std::streampos pos;

  bool secondev=false;
  bool bsp=false;
  JobType tjob=singlepoint;
  while ( std::getline( *m_file,line ) )
  {
    pos=m_file->tellg();
    if(line.find("DERIV_LEV")!=std::string::npos)
    {
      StringTokenizer token(line," \t");
      if( token.at(2) == "SECOND" ) secondev=true;
    }
#ifdef __GNUC__
#warning Is this really true
#endif
    if( line.find("GEO_METHOD") != std::string::npos)
    {

      bsp= ( line.find( "SINGLE_POINT" ) != std::string::npos );
      if(!bsp)
        tjob=opt;
      else
        if(secondev==true)
            tjob=freq;

      m_jobpos.push_back(JobHeader(tjob,pos));

    }

  }

  return m_jobpos;
}

