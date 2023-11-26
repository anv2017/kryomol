/*****************************************************************************************
                            pcmodelparser.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef PCMODELPARSER_H
#define PCMODELPARSER_H

#include "parser.h"
#include <map>

namespace kryomol
{
/**
A class for parsering of PCModel files
*/
class KRYOMOLPARSERS_API PCModelParser : public kryomol::Parser
{
    typedef std::map<const std::string, const int> AtomTable;
    typedef std::pair<const std::string, const int> AtomPair;

public:
    enum TypeTable { UNKNOWN, AMBER, MM2TEST, MM3, MMFF94, MMXCONST, OPLSAA};
    PCModelParser( const char* file);
    PCModelParser(std::istream* stream);
    ~PCModelParser();
    bool ParseFile(std::streampos pos=0);
private:
    void GetMolecule();
    bool GetConformer();
    void BuildAtomTable(TypeTable type);
    const AtomTable GetAtomTable();
    int ZFromString(const std::string symbol);

private:
  AtomTable m_atomtable;
};

}

#endif // PCMODELPARSER_H
