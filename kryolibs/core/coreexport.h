/*****************************************************************************************
                            coreexport.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef BUILDKRYOMOLCORE_H
#define BUILDKRYOMOLCORE_H


#if defined(WIN32) && !defined(KRYOMOLSTATIC)
#ifdef KRYOMOLCOREDLL
#define KRYOMOLCORE_EXPORT __declspec(dllexport)
#define KRYOMOLCORE_API __declspec(dllexport)
#else
#define KRYOMOLCORE_API __declspec(dllimport)
#endif
#else
#define KRYOMOLCORE_EXPORT
#define KRYOMOLCORE_API
#endif


#endif
