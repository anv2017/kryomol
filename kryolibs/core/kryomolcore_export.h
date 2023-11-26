/*****************************************************************************************
                            kryomolcore_export.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

/* needed for KDE_EXPORT and KDE_IMPORT macros */
#ifndef KRYOMOLCORE_EXPORT_H
#define KRYOMOLCORE_EXPORT_H

// needed for KDE_EXPORT and KDE_IMPORT macros
#ifdef KDEBUILD
#include <kdemacros.h>

#ifndef KRYOMOLCORE_EXPORT
# if defined(MAKE_KRYOMOLCORE_LIB)
   // We are building this library
#  define KRYOMOLCORE_EXPORT KDE_EXPORT
# else
   // We are using this library
#  define KRYOMOLCORE_EXPORT KDE_IMPORT
# endif
#endif

# ifndef KRYOMOLCORE_EXPORT_DEPRECATED
#  define KRYOMOLCORE_EXPORT_DEPRECATED KDE_DEPRECATED KRYOMOLCORE_EXPORT
# endif
#else
#define KRYOMOLCORE_EXPORT
#endif

#endif

