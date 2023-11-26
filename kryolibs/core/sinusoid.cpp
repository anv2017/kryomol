/*****************************************************************************************
                            sinusoid.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "sinusoid.h"


Sinusoid::Sinusoid()
{}

Sinusoid::~Sinusoid()
{}

//Allocates the data array of individual waves
bool Sinusoid::InitData(int n)
{
    m_fdata.resize(n);

    return true;
}



bool Sinusoid::NewData(int n)
{
    return true;
}

Sinusoid const &Sinusoid::operator=(  Sinusoid const & other)
{
    if (this!=&other)
    {
            m_flamplitude=other.m_flamplitude;
            m_zamplitude=other.m_zamplitude;
            m_fldecay=other.m_fldecay;
            m_flfrequency=other.m_flfrequency;
            m_flphase=other.m_flphase;
            m_flt1=other.m_flt1;
            m_J=other.m_J;
            m_multiplicity=other.m_multiplicity;
    }
    return(*this);
}

