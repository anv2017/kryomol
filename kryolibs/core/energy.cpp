/*****************************************************************************************
                            energy.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "energy.h"

using namespace kryomol;

constexpr double NA= 6.02214076e23;
constexpr double KB=3.166790920045901e-06;

Energy::unit Energy::Units() const { return m_units; }

void Energy::SetUnits(Energy::unit u) { m_units=u; }



double Energy::Value(unit u)
{
    if ( this->Units() == u ) return Value();
    return this->Value()*Factor(this->Units(),u);
}

double Energy::Kb(unit u /*=Hartree*/)
{
    double k=KB;
    if (u != HARTREE)
    {
        k*=Factor(HARTREE,u);
    }

    return k;
}
double Energy::Factor(Energy::unit u1, Energy::unit u2)
{
    if ( u1 == u2) return 1;
    bool invert =  u1 > u2;
    if (invert )
    {
        //std::swap is not compiling in the mac with clang 64
        unit tmp=u1;
        u1=u2;
        u2=tmp;
    }
    //std::swap(u1,u2);
    double factor;
    if ( u1 == HARTREE && u2 == EV )
        factor=27.2114;
    if ( u1 == HARTREE && u2 == J )
        factor =4.3597e-18*NA;
    if ( u1 == HARTREE && u2 == KJ )
        factor =4.3597e-21*NA;
    if ( u1 == HARTREE && u2 == CAL )
        factor=1.041993308181644261e-18*NA;
    if ( u1 == HARTREE && u2 == KCAL )
        factor=1.041993308181644261e-21*NA;
    if ( u1 == EV && u2 == J )
        factor=1.602160333349288422e-19*NA;
    if ( u1 == EV && u2 == KJ )
        factor=1.602160333349288422e-22*NA;
    if ( u1 == EV && u2 == CAL)
        factor=3.829255098827170949e-20*NA;
    if ( u1 == EV && u2 == KCAL )
        factor=3.829255098827170949e-23*NA;
    if ( u1 == J && u2 == KJ )
        factor = 1e-3;
    if ( u1 == J && u2 == CAL )
        factor=0.23900331477056468987;
    if ( u1 == J && u2 == KCAL)
        factor=0.23900331477056468987*1e-3;
    if ( u1 == KJ && u2 == CAL )
        factor=239.00331477056471385;
    if ( u1 == KJ && u2 == KCAL)
        factor=0.23900331477056471763;
    if ( u1 == CAL && u2 == KCAL)
        factor=1e-3;

    if ( invert ) factor=(1/factor);
    return factor;

}
