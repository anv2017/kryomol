/*****************************************************************************************
                            maestroparser.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef MAESTROPARSER_H
#define MAESTROPARSER_H

#include "parser.h"
#include "coordinate.h"


namespace kryomol
{

class KRYOMOLPARSERS_API MaestroParser : public kryomol::Parser
{
public:
  MaestroParser(const char* file);
  MaestroParser(std::istream* );
  ~MaestroParser();

  bool ParseFile(std::streampos pos=0);
private:
  void GetGeometries();
  void GetBlock(bool full);
  int AtomType(int type);
  std::vector< std::pair<size_t,int> > m_uatoms;
  std::vector<kryomol::Coordinate> m_newh;
  void AddHToCR3(size_t atom);
  void AddHToCR2(size_t atom);
  void AddHToCR(size_t atom);
  void AddHToCsp(size_t atom);
  void AddH2ToCsp2(size_t atom);
  void AddHToCsp2(size_t atom);
  void TransformUnitedAtoms();
  void TransformAtom(const std::pair<size_t,int>& atom);

};

}
#endif
