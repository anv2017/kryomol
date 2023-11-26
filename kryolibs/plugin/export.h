/*****************************************************************************************
                            export.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/


#ifndef PLUGINEXPORT_H
#define PLUGINEXPORT_H

#if defined(WIN32) && !defined(KRYOMOLSTATIC)
#include <windows.h>
#define KRYOMOL_EXPORT __declspec(dllexport)
#ifdef KRYOMOLDLL
#define KRYOMOL_API __declspec(dllexport)
#else
#define KRYOMOL_API __declspec(dllimport)
#endif
#else
#define KRYOMOL_API
#define KRYOMOL_EXPORT
#endif

#endif
