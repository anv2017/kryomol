/*****************************************************************************************
                            pdbparser.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <algorithm>
#include "pdbparser.h"
#include "molecule.h"
#include "stringtools.h"

using namespace kryomol;

PDBParser::PDBParser(const char* file)
: Parser(file)
{}



PDBParser::PDBParser(std::istream* stream) : Parser(stream)
{}

PDBParser::~PDBParser()
{}

bool PDBParser::ParseFile(std::streampos pos)
{
	
	m_file->seekg(0,std::ios::beg);
	std::string line;
	std::vector<std::streampos> spos;
	
	while(std::getline(*m_file,line))
	{
		if( line.find("MODEL") == 0 )
		{
			spos.push_back(m_file->tellg());
		}
	}
	if( spos.empty() )  //only one structure
	{ spos.push_back(0); }
	
	//Reserve the number of frames
	Molecules()->push_back(Molecule());
	
	Molecule& molecule=Molecules()->back();
	
	std::vector<std::streampos>::iterator it;
	
	
	it=spos.begin();
	m_file->clear();
        m_file->seekg(*it,std::ios::beg);
	while(std::getline(*m_file,line))
	{
		if ( line.find("HETATM")  != std::string::npos || line.find("ATOM") != std::string::npos)
		{
			std::vector<std::string> tokens=ParseAtomLine(line);
			
			if ( tokens.at(0) == "HETATM" || tokens.at(0) == "ATOM  " )
			{
                                //get atom symbol
				std::string symbol;
				
				//if we dont have symbol record extract it from the pdb atom name 
				if ( tokens.size() == 9 )
					symbol=tokens.at(8);
				else
					symbol=ExtractAtomName(tokens.at(2));
				
				Atom atom(symbol);
				atom.SetPDBName(tokens.at(2).c_str());

				
				std::vector<PDBResidue*>::iterator rt=molecule.Residues().begin();
				for(;rt!=molecule.Residues().end();++rt)
				{
					if ( (*rt)->Name() == tokens.at(3) && (*rt)->Index() == tokens.at(4) )
						break;
				}
				if ( rt != molecule.Residues().end() )
					atom.SetPDBResidue(*rt);
				else
				{
					PDBResidue* nresidue= new PDBResidue(tokens.at(3),tokens.at(4));
					molecule.Residues().push_back(nresidue);
					atom.SetPDBResidue(nresidue);
				}
				
				
				
				molecule.Atoms().push_back(atom);
			}
		}
		if( line.find("MODEL") == 0 ) break;
	}
	
	molecule.Frames().reserve(spos.size()); 

	for(it=spos.begin();it!=spos.end();++it)
        {
		m_file->clear();
		m_file->seekg(*it,std::ios::beg);
		//First count the atoms
		
		molecule.Frames().push_back(Frame(&molecule));
		Frame& frame=molecule.Frames().back();
		
		frame.XYZ().resize(molecule.Atoms().size());
		
		std::vector<Coordinate>::iterator ct=frame.XYZ().begin();
		
		while(std::getline(*m_file,line))
		{
#ifdef __GNUC__
#warning optimize this function call, redundant parsering
#endif
			if ( line.find("HETATM") != std::string::npos  || line.find("ATOM") != std::string::npos )
			{
				std::vector<std::string> tokens=ParseAtomLine(line);
				
				ct->x()=std::atof(tokens.at(5).c_str());
				ct->y()=std::atof(tokens.at(6).c_str());
				ct->z()=std::atof(tokens.at(7).c_str());
				++ct;
				
			}
			if( line.find("MODEL") == 0  ) break;
		}
		
		
	}
        return true;
	
}


std::vector<std::string> PDBParser::ParseAtomLine(const std::string& line)
{
	std::vector<std::string> tokens;
	tokens.reserve(9);
	tokens.push_back(line.substr(0,6)); //"ATOM  " or "HETATOM" 0
	tokens.push_back(line.substr(6,5)); //Serial number 1
	
	std::string pdbname=line.substr(12,4);
	
	
	pdbname.erase(std::remove(pdbname.begin(), pdbname.end(),' '), pdbname.end()); 
	tokens.push_back(pdbname);// atom pdb name 2
	
	std::string resname=line.substr(17,3);
	resname.erase(std::remove(resname.begin(),resname.end(),' '),resname.end());
	tokens.push_back(resname);// atom residue name
	
	
	std::string resindex=line.substr(22,5);
	resindex.erase(std::remove(resindex.begin(),resindex.end(),' '),resindex.end());
	tokens.push_back(resindex);//resnumber is not an index see especification!!! 4
	
	tokens.push_back(line.substr(30,8));//x coordinate 5
	tokens.push_back(line.substr(38,8));//y coordinate 6
	tokens.push_back(line.substr(46,8));//z coordinate 7
	
       //we may not have atomic number here
	if (line.size() >= 78 )
	{
		std::string symbol=line.substr(76,2);
		symbol.erase(std::remove(symbol.begin(),symbol.end(),' '),symbol.end());
		tokens.push_back(symbol);//element symbol 8
	}
	
	return tokens;
	
}

std::string PDBParser::ExtractAtomName(const std::string& s)
{
	PeriodicTable::const_iterator it;
	std::vector<PeriodicTable::const_iterator> atomnames;
	std::vector<size_t> cpos;
	for(it=GetPeriodicTable().begin();it!=GetPeriodicTable().end();++it)
	{
		size_t p=s.find(const_cast<std::string&>(it->first),0);
                if ( p != std::string::npos )
		{
			atomnames.push_back(it);
			cpos.push_back(s.find(p,0));
		}
		
		if( atomnames.size() == 1 ) break;
	}
	
	if ( !atomnames.empty() ) 
		return atomnames.at(0)->first;
	else return "X";
	
}
