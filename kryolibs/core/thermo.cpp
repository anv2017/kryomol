/*****************************************************************************************
                            thermo.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <limits>

#include "thermo.h"
#include "molecule.h"

using namespace kryomol;

Thermostat::Thermostat() :  m_ensamble(NVE), m_temperature(298.15)
{}


Thermostat::~Thermostat()
{}

void Thermostat::SetPopulations(Molecule& molecule) const
{
  switch( m_ensamble)
  {
  case NVT:  //Boltzmann populations
      std::cout << "Setting boltzmann populations" << std::endl;
    SetBoltzmannPopulations(molecule);
    break;
  default:
    std::cout << "Setting populations" << std::endl;
    molecule.Populations()=std::vector<double>(molecule.Frames().size(),1/static_cast<double>(molecule.Frames().size()));
    break;
  }
}

void Thermostat::SetBoltzmannPopulations(Molecule& molecule) const

{
  double k=1.98588e-3; //kcal/mol

  std::vector<double> pop(molecule.Frames().size());
  std::vector<double>::iterator it;
  std::vector<Frame>::const_iterator ft=molecule.Frames().begin();

  for(;ft!=molecule.Frames().end();++ft)
  {
      if ( !ft->PotentialEnergy() )
      {
          throw kryomol::Exception("Potential energy not defined");
      }
  }

  ft=molecule.Frames().begin();
  double z=0;
  for(it=pop.begin();it!=pop.end();++it,++ft)
  {

    double e= ft->PotentialEnergy().Value();
    (*it)=exp( -e/(k*m_temperature));
    z+=(*it);
  }

  for(it=pop.begin();it!=pop.end();++it)
  {
    (*it)/=z;
  }
  molecule.Populations()=pop;

}
