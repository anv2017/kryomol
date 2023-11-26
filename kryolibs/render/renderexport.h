/*****************************************************************************************
                            renderexport.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef BUILDQRYOMOLRENDER_H
#define BUILDQRYOMOLRENDER_H


#if defined(WIN32) && !defined(KRYOMOLSTATIC)
#ifdef KRYOMOLRENDERDLL
#define KRYOMOLRENDER_EXPORT __declspec(dllexport)
#define KRYOMOLRENDER_API __declspec(dllexport)
#else
#define KRYOMOLRENDER_API __declspec(dllimport)
#endif
#else
#define KRYOMOLRENDER_EXPORT
#define KRYOMOLRENDER_API
#endif

#endif // RENDEREXPORT_H
