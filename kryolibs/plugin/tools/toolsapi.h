#ifndef TOOLSAPI_H
#define TOOLSAPI_H

#ifdef WIN32
#ifdef BUILDTOOLSDLL
#define TOOLS_API __declspec(dllexport)
#else
#define TOOLS_API __declspec(dllimport)
#endif
#else
#define TOOLS_API 
#endif

#endif
