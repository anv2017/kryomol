/*****************************************************************************************
                            utilities.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef UTILITIES_H
#define UTILITIES_H
#include <math.h>
inline float MulDiv(float a, float b,float c){
  return floor((a*b/c)+0.5f);
}

#endif
