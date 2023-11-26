/*****************************************************************************************
                            quantumcoupling.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QUANTUMCOUPLING_H
#define QUANTUMCOUPLING_H


#include "coreexport.h"
#include "mathtools.h"

namespace kryomol
{ 
  class QuantumCouplingPrivate;

  class KRYOMOLCORE_API QuantumCoupling
  {
    public:
    QuantumCoupling(size_t i, size_t j);
    QuantumCoupling(size_t i,size_t j, double v);
    QuantumCoupling( const QuantumCoupling& c );
    QuantumCoupling& operator= ( const QuantumCoupling& c );
    ~QuantumCoupling();
    const size_t& I() const;
    const size_t& J() const;
    size_t& I();
    size_t& J();
    const Scalar<double>& Value() const;
    void SetValue(double v);
    const Scalar<double>& FermiContact() const;
    const Scalar<double>& SpinDipole() const;
    const Scalar<double>& FCSDCrossTerm() const;
    const Scalar<double>& ParamagneticSpinOrbit() const;
    const Scalar<double>& DiamagneticSpinOrbit() const;
    void SetFermiContact(double v);
    void SetSpinDipole(double v);
    void SetFCSDCrossTerm(double v);
    void SetParamagneticSpinOrbit(double v);
    void SetDiamagneticSpinOrbit(double v);
    private:
    QuantumCouplingPrivate* _d;
  };
}
#endif //QUANTUMCOUPLING_H
