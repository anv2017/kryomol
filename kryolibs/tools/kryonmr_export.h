/* needed for KDE_EXPORT and KDE_IMPORT macros */
#ifndef KRYOMOL_EXPORT_H
#define KRYOMOL_EXPORT_H

// needed for KDE_EXPORT and KDE_IMPORT macros
#include <kdemacros.h>

#ifndef KRYOMOL_EXPORT
# if defined(MAKE_KRYOMOL_LIB)
   // We are building this library
#  define KRYOMOL_EXPORT KDE_EXPORT
# else
   // We are using this library
#  define KRYOMOL_EXPORT KDE_IMPORT
# endif
#endif

# ifndef KRYOMOL_EXPORT_DEPRECATED
#  define KRYOMOL_EXPORT_DEPRECATED KDE_DEPRECATED KRYOMOL_EXPORT
# endif

#endif

