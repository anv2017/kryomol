/*****************************************************************************************
                            chemicalshift.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "chemicalshift.h"

bool kryomol::operator==(const kryomol::ChemicalShift & a , const  kryomol::ChemicalShift & b)
{
  return (  a.I() == b.I() );
}

using namespace kryomol;

ChemicalShift::ChemicalShift(size_t i) : _first(i), _hasexperimentalvalue(false)
{
}

ChemicalShift::ChemicalShift(size_t i,  double ev) : _first(i),  _experimentalvalue(ev), _hasexperimentalvalue(true)
{
}
ChemicalShift::ChemicalShift(size_t i,double ev,double err) : _first(i), _experimentalvalue(ev), _hasexperimentalvalue(true)
{
  SetStandardDeviation(err);
}

 void ChemicalShift::SetExperimentalValue(double ev) { _experimentalvalue=ev; _hasexperimentalvalue=true; 
}





