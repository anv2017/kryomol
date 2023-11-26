/*****************************************************************************************
                            thermo.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QRYOMOLTHERMO_H
#define QRYOMOLTHERMO_H

#include <vector>
#include "coreexport.h"

/** */
namespace kryomol {

class Molecule;

KRYOMOLCORE_API enum ensamble { NVE=0, NVT };

/**
@brief A class to manage thermodynamics

*/
class KRYOMOLCORE_API Thermostat {
public:
    Thermostat();
    virtual ~Thermostat ();
    /** Set the temperature to t Kelvin*/
    virtual void SetTemperature(double t) { m_temperature=t; }
    /** @return the temperature T of the thermostat*/
    double Temperature() const { return m_temperature; }
//    void SetMolecule(Molecule* mol) { m_molecule=mol; }
    /** set the thermodynamic ensamble to e , either energy constant (NVE) or temperature constant (NVT*)*/
    virtual void SetEnsamble(ensamble e) { m_ensamble=e; }
    /** @return the  kind of thermodynamic ensamble , either energy constant (NVE) or temperature constant (NVT*)*/
    ensamble Ensamble() const { return m_ensamble; }
    /** @return populations of all conformers for the @see qryomol::Molecule molecule*/
    virtual void SetPopulations(Molecule& molecule) const;
private: //private methodsAcesParser
   void SetBoltzmannPopulations(Molecule& molecule) const;
private:
  ensamble m_ensamble;
  double m_temperature;
};



}

#endif
