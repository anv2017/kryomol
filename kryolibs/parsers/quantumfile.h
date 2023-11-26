/***************************************************************************
                          quantumfile.h  -  description
                             -------------------
    copyright            : (C) 2007 by A. Navarro-Vazquez
    email                : qoajnv@usc.es
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef QUANTUMFILE_H
#define QUANTUMFILE_H

/**
@author Armando Navarro-Vazquez
*/
#include <fstream>
#include <istream>

/**This class try to guess the ouput file type */
class QuantumFile{
public:
	QuantumFile(const char* file) { m_file = new std::ifstream(file,std::ios::binary); m_bstreamcreated=true; }
    QuantumFile(std::istream* stream) : m_file(stream) { m_bstreamcreated= false; }

    ~QuantumFile() { if(m_bstreamcreated) delete m_file; };
    
    enum qfiletype { None, Gaussian, GaussianArchive, GaussianInput, Aces, Gamess, Mol, Psi3, PDB, MacroModel,Maestro,XYZ };
    
    qfiletype GetFileType();
    
private:
  std::istream* m_file;
  bool m_bstreamcreated;

};

#endif
