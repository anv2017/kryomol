/*****************************************************************************************
                            parserfactory.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef PARSERFACTORY_H
#define PARSERFACTORY_H

#include "parsersexport.h"

#include <fstream>
#include <istream>

#ifdef __MINGW32__
#include <filesystem>
#endif

namespace kryomol
{
class MagicFile;
class Parser;

/**A class to build specialized parsers */
class KRYOMOLPARSERS_API ParserFactory
{
public:
    ParserFactory ( const char* file );
    ParserFactory ( std::istream* stream );
#ifdef __MINGW32__
    ParserFactory( std::filesystem::path p);
#endif
    ~ParserFactory();
    /** return a pointer to a specialized parser, NULL if file type cannot be discerned*/
    Parser* BuildParser();
    enum filetype { None, Aces, Gaussian, GaussianArchive, GaussianInput, GaussianFile, GaussianCube, Gamess, MdlV2000, PDB, MacroModel, Maestro, NwChem, XYZ, HyperChem, PCModel, Orca };
    filetype GetFileType();
    bool isGaussianFile();
    bool existDensity();
    bool existOrbitals();
    bool existAlphaBeta();
private:
    std::istream* m_stream;
    bool m_bstreamcreated;

};

}
#endif
