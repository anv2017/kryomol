/*****************************************************************************************
                            sinusoid.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#if !defined(AFX_SINUSOID_H__48A6CB59_D7BA_4FD2_BE7F_57F6CC68CC16__INCLUDED_)
#define AFX_SINUSOID_H__48A6CB59_D7BA_4FD2_BE7F_57F6CC68CC16__INCLUDED_

#include "fidarray.h"

class Sinusoid
{
public:

    Sinusoid();
    ~Sinusoid();

public:
    bool NewData(int n);
    int m_multiplicity;
    double m_zamplitude;
    Sinusoid const & operator=(Sinusoid const  & other);
    bool InitData(int n);
    float m_flphase;  //phase of individual fid
    float m_flamplitude; //amplitude of individual fid
    float m_fldecay;  //decay time of individual fid
    float m_flfrequency; //frequency of individual fid
    float m_flt1;   //T1 of each frequency
    float m_J; //coupling constant in Hz
    fidarray m_fdata; //array of points

};

#endif // !defined(AFX_SINUSOID_H__48A6CB59_D7BA_4FD2_BE7F_57F6CC68CC16__INCLUDED_)
