/***************************************************************************
                          quantumfile.cpp  -  description
                             -------------------
    copyright            : (C) 2007 by A. Navarro-Vazquez
    email                : qoajnv@usc.es
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "quantumfile.h"
#include <iostream>
#include "stringtools.h"

using namespace nmrdev;
QuantumFile::qfiletype QuantumFile::GetFileType()
{

  //first try gaussian
  std::string line;

  while(std::getline(*m_file,line))
  {
    if(line.find("Standard orientation",0) != std::string::npos ||  line.find("Input orientation",0) != std::string::npos )
    {
     std::cout << "Gaussian Output File" << std::endl;
      return Gaussian;
     }
  }

  m_file->clear();
  m_file->seekg(0,std::ios::beg);
  std::getline(*m_file,line);
  StringTokenizer tok(line," \t\r,");
  if( tok.size()== 2 )
  {
    if ( isinteger(tok.at(0)) && isinteger(tok.at(1)) )
      return GaussianInput;
  }

  m_file->clear();
  m_file->seekg(0,std::ios::beg);
  while(std::getline(*m_file,line))
  {
    if(line.find("Z-matrix") != std::string::npos)
      return Aces;
  }

  m_file->clear();
  m_file->seekg(0,std::ios::beg);
  while(std::getline(*m_file,line))
  {
    if(line.find("GAMESS VERSION") != std::string::npos)
      return Gamess;
  }

  m_file->clear();
  m_file->seekg(0,std::ios::beg);
  while(std::getline(*m_file,line))
  {
    if(line.find("s_m_m2io_version") != std::string::npos)
    {
      std::cout << "MacroModel Maestro format" << std::endl;
      return Maestro;
    }

  }

  m_file->clear();
  m_file->seekg(0,std::ios::beg);

  //well, I dont know exactly what to search for, maybe 5 integers?
  while(std::getline(*m_file,line) )
  {
    bool isheader=false;
    StringTokenizer token(line," \t");
    if ( token.size() > 4 && token.size() < 8)
    {
      isheader=true;
      for(size_t i=0;i<5;i++)
      {
        bool kk= isinteger(token.at(i));
        //std::cout << "token is " << token.at(i) << "res=" << kk << std::endl;
        if( !kk) isheader=false;
      }
    }

    if( isheader ) return Mol;

  }


  m_file->clear();
  m_file->seekg(0,std::ios::beg);
  while(std::getline(*m_file,line))
  {
    if(line.find("WELCOME") != std::string::npos)
    {
      std::getline(*m_file,line);
      if(line.find("PSI") != std::string::npos) return Psi3;
    }
  }

  m_file->clear();
  m_file->seekg(0,std::ios::beg);
  std::getline(*m_file,line);
  if( ( line.find( "1\\1\\" ) != std::string::npos  ) ||  ( line.find( "1|1|") != std::string::npos ) ) return GaussianArchive;

  m_file->clear();
  m_file->seekg(0,std::ios::beg);
  while(std::getline(*m_file,line) )
  {
    if ( (line.find( "HETATM" ) != std::string::npos ) || line.find( "ATOM" ) != std::string::npos )
      return PDB;
  }

  //   //MacroModel dynamics
  //   m_file->clear();
  //   m_file->seekg(0,std::ios::beg);
  //   std::getline(*m_file,line);
  //   StringTokenizer tok1(line," \t");
  //   std::getline(*m_file,line);
  //   StringTokenizer tok2(line," \t");
  //   if( tok1.size() == 5 && tok2.size() > 10 )
  //    return MacroModelDyn;
  //MacroModel dynamics
  m_file->clear();
  m_file->seekg(0,std::ios::beg);
  std::getline(*m_file,line);
  std::getline(*m_file,line);
  {
    StringTokenizer token(line," \t\r");
    if (token.size() > 13 )
    {
      for(int i=0;i<13;i++)
      {
        if ( !isinteger(token.at(i))) goto lb;
      }
      std::cout << "MacroModel Old Format" << std::endl;
      return MacroModel;
    }
  }
lb:
  m_file->clear();
  m_file->seekg(0,std::ios::beg);
  std::getline(*m_file,line);
  StringTokenizer token(line," \t\r");
  if( token.size() == 1 )
    if( nmrdev::isnum(token.at(0)) ) return XYZ;

  return None;

}



