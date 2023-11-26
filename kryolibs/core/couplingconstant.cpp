/*****************************************************************************************
                            couplingconstant.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <iostream>

#include "couplingconstant.h"


bool kryomol::operator==(const CouplingConstant & a , const  CouplingConstant & b)
{
  return ( ( a.I() == b.I() ) && ( a.J() == b.J() ) );
}
std::ostream& kryomol::operator << ( std::ostream& s,const CouplingConstant& c )
{
  s << "(" << (c.I()+1) << "," << (c.J()+1) << ")=" << c.ExperimentalValue() << std::endl;
  return s;
}
using namespace kryomol;

CouplingConstant::CouplingConstant(size_t i, size_t j) : _first(i), _second(j), _hasexperimentalvalue(false)
{
}

CouplingConstant::CouplingConstant(size_t i, size_t j, double ev) : _first(i), _second(j), _experimentalvalue(ev), _hasexperimentalvalue(true)
{
}
CouplingConstant::CouplingConstant(size_t i, size_t j, double ev,double err) : _first(i), _second(j), _experimentalvalue(ev), _hasexperimentalvalue(true)
{
  SetStandardDeviation(err);
}

 void CouplingConstant::SetExperimentalValue(double ev) { _experimentalvalue=ev; _hasexperimentalvalue=true; 
}

