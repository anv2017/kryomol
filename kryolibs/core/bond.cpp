/*****************************************************************************************
                            bond.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/
 
#include "bond.h"

using namespace kryomol;

Bond::Bond ( size_t i,size_t j, order e ) : m_order ( e )
{
  if ( i < j )
  {
    m_i=i;
    m_j=j;
  }
  else
  {
    m_j=i;
    m_i=j;
  }
}
