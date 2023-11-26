/*****************************************************************************************
                            mathtools.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <string>

#include "mathtools.h"

void kryomol::MathInfo(std::string& library, std::string& storage)
{

#ifdef WITH_MKL
		library="Intel MKL";
		storage="column_major";
		return;
#endif
	library="undefined";
	storage="undefined";
	return;
}



