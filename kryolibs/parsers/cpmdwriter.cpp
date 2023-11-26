/*****************************************************************************************
                            cpmdwriter.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <iomanip>
#include "cpmdwriter.h"
#include "molecule.h"
#include <iostream>
#include <vector>

using namespace kryomol;

CPMDWriter::CPMDWriter(Molecule* molecule): m_molecule (molecule)
{}

CPMDWriter::~CPMDWriter()
{}

void kryomol::operator << ( std::ostream& s, const CPMDWriter& w )
{
    std::cout << "Exporting CPMD input file" << std::endl;
 s << "&CPMD" << std::endl << "&END" << std::endl;
 s << "&SYSTEM" << std::endl;
 s << "ANGSTROM" << std::endl;
 s << "&END" << std::endl;
 s << "&ATOMS" << std::endl;
 //OK lets do groups
 std::vector<Atom>::const_iterator at;
 std::vector <std::string> names;
 std::vector<std::string>::const_iterator st;
 for ( at=w.m_molecule->Atoms().begin();at!=w.m_molecule->Atoms().end();++at )
 {

   for ( st=names.begin();st!=names.end();++st )
   {
     if ( ( *st ) == at->Symbol() ) break;
   }
   if ( ! ( st != names.end() ) ) names.push_back ( at->Symbol() );
 }
 for ( st=names.begin();st!=names.end();++st )
 {
   s << "*" << ( *st ) << ".psp" << std::endl;
   s << "LMAX=" << std::endl;
   size_t counter=0;
   for ( at=w.m_molecule->Atoms().begin();at!=w.m_molecule->Atoms().end();++at )
   {
     if ( at->Symbol() == ( *st ) ) ++counter;
   }
   s << " " << counter << std::endl;

   std::vector<Coordinate>::const_iterator ft = w.m_molecule->CurrentFrame().XYZ().begin();
   for ( at=w.m_molecule->Atoms().begin();at!=w.m_molecule->Atoms().end();++at,++ft )
   {
     if ( at->Symbol() == ( *st ) )
     {
       s << std::setw ( 10 ) << std::setiosflags ( std::ios::right )
       << std::setiosflags ( std::ios::fixed ) << std::setprecision ( 5 )
         << ft->x() << std::setw ( 10 ) << ft->y() << std::setw ( 10 ) << ft->z()
         << std::resetiosflags ( std::ios::right ) << std::endl;
     }
   }
 }
 s << "&END" << std::endl;
}

