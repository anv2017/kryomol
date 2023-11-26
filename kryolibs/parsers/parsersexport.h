/*****************************************************************************************
                            parsersexport.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef BUILDKRYOMOLPARSERS_H
#define BUILDKRYOMOLPARSERS_H


#if defined(WIN32) && !defined(KRYOMOLSTATIC)
#ifdef KRYOMOLPARSERS_DLL
#define KRYOMOLPARSERS_EXPORT __declspec(dllexport)
#define KRYOMOLPARSERS_API __declspec(dllexport)
#else
#define KRYOMOLPARSERS_API __declspec(dllimport)
#endif
#else
#define KRYOMOLPARSERS_EXPORT
#define KRYOMOLPARSERS_API
#endif


#endif
