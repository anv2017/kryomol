/*****************************************************************************************
                            parser.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <clocale>
#include "parser.h"


using namespace kryomol;

Parser::Parser ( const char* inputfile ) : m_bcreated ( true )
{
    m_file= new std::ifstream ( inputfile );
}

Parser::Parser ( std::istream* stream ) : m_file ( stream ) , m_bcreated ( false )
{
}

Parser::~Parser()
{
    if ( m_bcreated ) delete m_file;
}


void Parser::ParseChemicalShifts()
{
}
void Parser::ParseCouplingConstants(std::vector<QuantumCoupling>& c)
{
}

D2Array<double>  Parser::ParseMagneticSusceptibility()
{
    return D2Array<double>();
}

/** Parse will set the locale to C, call the actual virtual implementation ParseFile and restore the locale*/
void Parser::Parse(std::streampos pos)
{
    const char* locale=std::setlocale(LC_NUMERIC,"C");
    ParseFile(pos);
    std::setlocale(LC_NUMERIC,locale);
}

std::vector<JobHeader>& Parser::Jobs()
{
  if ( m_jobpos.empty() ) m_jobpos.push_back(JobHeader(singlepoint,0));
  return m_jobpos;
}
