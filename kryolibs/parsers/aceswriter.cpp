/*****************************************************************************************
                            aceswriter.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "aceswriter.h"
#include "molecule.h"
#include <iostream>
#include <vector>
#include <iomanip>
#include "stringtools.h"

using namespace kryomol;

AcesWriter::AcesWriter(Molecule* molecule): m_molecule (molecule)
{}

AcesWriter::~AcesWriter()
{
}

std::ostream& kryomol::operator << (std::ostream& s, const AcesWriter& aw)
{
    std::cout << "Exporting ACES input file" << std::endl;

   Molecule* molecule=aw.m_molecule;

   std::stringstream svalues;
   size_t counter;
   size_t at;

   for(at=0, counter=0; at<molecule->Atoms().size();at++,counter++)
   {

       s << toupper(molecule->Atoms().at(at).Symbol()) << " ";

        if( counter )
        {
         AcesWriter::zmatlabel r;
         r.label << "r" <<counter;
         r.value=Coordinate::Distance(molecule->CurrentFrame().XYZ().at(at),molecule->CurrentFrame().XYZ().at(at-1));
         s <<  counter << " " << r.label.str()  << " ";

         svalues <<  r.label.str() << "=" <<  std::setiosflags(std::ios::fixed) <<std::setprecision(4) <<
           r.value   << std::resetiosflags(std::ios::fixed) << std::endl;
         }

         if( counter > 1)
         {
           AcesWriter::zmatlabel a;
            a.label << "a" << counter-1;
            a.value=Coordinate::GetAngle(molecule->CurrentFrame().XYZ().at(counter),molecule->CurrentFrame().XYZ().at(counter-1),molecule->CurrentFrame().XYZ().at(counter-2),true);
            s << counter-1 << " " << a.label.str() << " ";
            svalues << a.label.str() << "=" <<  std::setiosflags(std::ios::fixed) <<std::setprecision(3) <<
           a.value   << std::resetiosflags(std::ios::fixed) << std::endl;
            }

         if(counter > 2)
         {
           AcesWriter::zmatlabel d;
           d.label << "d" << counter-2;
           d.value=Coordinate::GetDihedral(molecule->CurrentFrame().XYZ().at(counter),molecule->CurrentFrame().XYZ().at(counter-1),molecule->CurrentFrame().XYZ().at(counter-2),molecule->CurrentFrame().XYZ().at(counter-3),true);
           s << counter -2 << " " << d.label.str();
           svalues << d.label.str() << "=" << std::setiosflags(std::ios::fixed) <<std::setprecision(2) <<
           d.value   << std::resetiosflags(std::ios::fixed) << std::endl;
         }

         s << std::endl;


   }

   s << std::endl << svalues.str() << std::endl;

   return s;


}
