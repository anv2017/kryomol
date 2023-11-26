/* needed for KDE_EXPORT and KDE_IMPORT macros */
#ifndef KRYOTOOLS_EXPORT_H
#define KRYOTOOLS_EXPORT_H

// needed for KDE_EXPORT and KDE_IMPORT macros
#include <kdemacros.h>

#ifndef KRYOTOOLS_EXPORT
# if defined(MAKE_KRYOTOOLS_LIB)
   // We are building this library
#  define KRYOTOOLS_EXPORT KDE_EXPORT
# else
   // We are using this library
#  define KRYOTOOLS_EXPORT KDE_IMPORT
# endif
#endif

# ifndef KRYOTOOLS_EXPORT_DEPRECATED
#  define KRYOTOOLS_EXPORT_DEPRECATED KDE_DEPRECATED KRYOTOOLS_EXPORT
# endif

#endif

