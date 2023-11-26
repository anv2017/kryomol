/*****************************************************************************************
                            url.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QRYOMOLURL_H
#define QRYOMOLURL_H

class QString;

#include "export.h"

namespace kryomol {

  KRYOMOL_API void OpenURL(const QString& url);

};

#endif
