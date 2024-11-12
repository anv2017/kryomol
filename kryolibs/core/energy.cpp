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

Energy::unity Energy::Unity() const { return m_unity; }

void Energy::SetUnity(Energy::unity u) { m_unity=u; }



double Energy::Value(unity u)
{
    if ( this->Unity() == u ) return Value();
    return this->Value()*Factor(this->Unity(),u);
}


double Energy::Factor(Energy::unity u1, Energy::unity u2)
{
    if ( u1 == u2) return 1;
    bool invert =  u1 > u2;
    if (invert )
    {
        //std::swap is not compiling in the mac with clang 64
        unity tmp=u1;
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
