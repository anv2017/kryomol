/*****************************************************************************************
                            kryoparsers_export.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

/* needed for KDE_EXPORT and KDE_IMPORT macros */
#ifndef KRYOPARSERS_EXPORT_H
#define KRYOPARSERS_EXPORT_H

// needed for KDE_EXPORT and KDE_IMPORT macros
#include <kdemacros.h>

#ifndef KRYOPARSERS_EXPORT
# if defined(MAKE_KRYOPARSERS_LIB)
   // We are building this library
#  define KRYOPARSERS_EXPORT KDE_EXPORT
# else
   // We are using this library
#  define KRYOPARSERS_EXPORT KDE_IMPORT
# endif
#endif

# ifndef KRYOPARSERS_EXPORT_DEPRECATED
#  define KRYOPARSERS_EXPORT_DEPRECATED KDE_DEPRECATED KRYOPARSERS_EXPORT
# endif

#endif

