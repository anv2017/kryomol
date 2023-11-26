/*****************************************************************************************
                            qdoubleslider.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "qdoubleslider.h"
#include "math.h"

void QDoubleSlider::OnValueChanged(int value)
{
    m_divisions = this->maximum()-this->minimum()+1;
    m_singlestep = (log(m_max)-log(m_min))/(m_divisions-1);
    m_value = 1/exp(value*m_singlestep);

    emit valueChanged(m_value);
}

void QDoubleSlider::OnSliderReleased()
{
    emit valueChanged(m_value);
}
