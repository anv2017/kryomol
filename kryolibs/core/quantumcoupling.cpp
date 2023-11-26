/*****************************************************************************************
                            quantumcoupling.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "quantumcoupling.h"

class kryomol::QuantumCouplingPrivate
{
    public:
    QuantumCouplingPrivate() {}
    ~QuantumCouplingPrivate() {}
    Scalar<double> totalvalue;
    Scalar<double> fermicontact;
    Scalar<double> spindipole;
    Scalar<double> crossfcsd;
    Scalar<double> diaspinorbit;
    Scalar<double> paraspinorbit;
    size_t I;
    size_t J;
};

using namespace kryomol;

QuantumCoupling::QuantumCoupling ( size_t i, size_t j )
{
  _d = new QuantumCouplingPrivate();
  if ( i < j )
  {
  _d->I=i;
  _d->J=j;
  }
  else
  {
    _d->I=j;
    _d->J=i;
  }
}

QuantumCoupling::QuantumCoupling ( size_t i, size_t j, double v )
{
  _d = new QuantumCouplingPrivate();
  if ( i < j )
  {
  _d->I=i;
  _d->J=j;
  }
  else
  {
    _d->I=j;
    _d->J=i;
  }
  _d->totalvalue=v;
}

QuantumCoupling::~QuantumCoupling()
{
  delete _d;
}

QuantumCoupling::QuantumCoupling ( const QuantumCoupling& c )
{
  _d = new QuantumCouplingPrivate();
  _d->I=c.I();
  _d->J=c.J();
  _d->totalvalue=c.Value();
  _d->fermicontact=c.FermiContact();
  _d->spindipole=c.SpinDipole();
  _d->crossfcsd=c.FCSDCrossTerm();
  _d->paraspinorbit=c.ParamagneticSpinOrbit();
  _d->diaspinorbit=c.DiamagneticSpinOrbit();

}

QuantumCoupling& QuantumCoupling::operator= ( const QuantumCoupling& c )
{
  if ( &c != this )
  {
    delete _d;
    _d = new QuantumCouplingPrivate();
    _d->I=c.I();
    _d->J=c.J();
    _d->totalvalue=c.Value();
    _d->fermicontact=c.FermiContact();
    _d->spindipole=c.SpinDipole();
    _d->crossfcsd=c.FCSDCrossTerm();
    _d->paraspinorbit=c.ParamagneticSpinOrbit();
    _d->diaspinorbit=c.DiamagneticSpinOrbit();
  }
  return *this;
}

const size_t& QuantumCoupling::I() const
{
  return _d->I;
}

const size_t& QuantumCoupling::J() const
{
  return _d->J;
}

size_t& QuantumCoupling::I()
{
  return _d->I;
}

size_t& QuantumCoupling::J() 
{
  return _d->J;
}

const Scalar<double>& QuantumCoupling::Value() const
{
  return _d->totalvalue;
}

void QuantumCoupling::SetValue(double v)
{
  _d->totalvalue=v;
}

const Scalar<double>& QuantumCoupling::FermiContact() const
{
  return _d->fermicontact;
}
    
const Scalar<double>& QuantumCoupling::SpinDipole() const
{
  return _d->spindipole;
}

const Scalar<double>& QuantumCoupling::FCSDCrossTerm() const
{
  return _d->crossfcsd;
}

const Scalar<double>& QuantumCoupling::ParamagneticSpinOrbit() const
{
  return _d->paraspinorbit;
}

const Scalar<double>& QuantumCoupling::DiamagneticSpinOrbit() const
{
  return _d->diaspinorbit;
}

void QuantumCoupling::SetFermiContact(double v)
{
  _d->fermicontact=v;
}

void QuantumCoupling::SetSpinDipole(double v)
{
  _d->spindipole=v;
}

void QuantumCoupling::SetFCSDCrossTerm(double v)
{
  _d->crossfcsd=v;
}

void QuantumCoupling::SetParamagneticSpinOrbit(double v)
{
  _d->paraspinorbit=v;
}

void QuantumCoupling::SetDiamagneticSpinOrbit(double v)
{
  _d->diaspinorbit=v;
}
