/***************************************************************************
                          pdbparser.cpp  -  description
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

#include "pdbparser.h"
#include "molecule.h"
#include "stringtools.h"
using namespace nmrdev;

PDBParser::PDBParser(const char* file)
    : Parser(file)
{}



PDBParser::PDBParser(std::istream* stream) : Parser(stream)
{}

PDBParser::~PDBParser()
{}

void PDBParser::ParseFile(std::streampos pos)
{


  m_file->seekg(0,std::ios::beg);
  std::string line;
  std::vector<std::streampos> spos;
  
  //spos.push_back(0);
  
  
  while(std::getline(*m_file,line))
  {
    if( line.find("MODEL") != std::string::npos)
    {
       //if ( std::getline(*m_file,line) )
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
 // bool bconected=false;
  while(std::getline(*m_file,line))
  {
    StringTokenizer token(line);
    size_t i;
    for(i=0;i<token.size();i++)
    {
      if(token[i]=="HETATM" || token[i]=="ATOM") //Do the stuff
      {
        int n=atoi(token.at(i+1).c_str())-1;
        std::string orname=token.at(i+2);
	std::string name=orname;
        std::string nn=ExtractAtomName(remove_numbers(name));
        std::string residue=token.at(i+3);
        int skipchain=0;
     
        if ( !nmrdev::isinteger(token.at(i+4) ) ) 
        {
           if (  token.at(i+4).size()  == 1 ) ++skipchain;
           else token.at(i+4)= token.at(i+4).substr(1);
        }
      
      //  std::cout << token.at(i+4+skipchain) << std::endl;

        int rnumber=atoi(token.at(i+4+skipchain).c_str());
        std::vector<PDBResidue*>::iterator rt;
        Atom atom(nn);
        atom.SetPDBName(orname);
        std::cout << name << " " << residue << " " << rnumber << std::endl;
        for(rt=molecule.Residues().begin();rt!=molecule.Residues().end();++rt)
        {
          if ( (*rt)->Name() == residue && (*rt)->Index() == rnumber )
             break;
        }
        if ( rt != molecule.Residues().end() )
        {
          atom.SetPDBResidue(*rt);
         
        }
        else
        {
          PDBResidue* nresidue= new PDBResidue(residue,rnumber);
          molecule.Residues().push_back(nresidue);
          atom.SetPDBResidue(nresidue);
        }

        molecule.Atoms().push_back(atom);
 
      }

     



      /*  if(token[i]=="CONECT")
        {
          bconected=true;
          int natom=atoi(token.at(i+1).c_str());
          int nconexions=0;
          std::vector<std::string> s;
          s.reserve(8);


          while((i+2+nconexions)<token.size())
          {
            if(token.at(i+2+nconexions)=="\r")
              break;
            s.push_back(token.at(i+2+nconexions));
            nconexions++;

          }


          for(int j=0;j<nconexions;j++)
          {

            molecule.m_atom[natom-1].conexion.push_back( bond( atoi(s.at(j).c_str() )-1) );
          }
          //Well, we should use the data from CONNECT but by now simply use our algorithm
          molecule.SetBondOrders();
        }




      }*/
    }
    if( line.find("MODEL") != std::string::npos /*|| line.find("END") != std::string::npos*/ ) break;
  }


   molecule.Frames().reserve(spos.size()); 

  int structure=0;
  for(it=spos.begin();it!=spos.end();++it)
  {
    //CMolecule molecule;
    std::cout << "Structure" <<  ++structure << std::endl;
    int atomnumber=0;
    m_file->clear();
    m_file->seekg(*it,std::ios::beg);
    //First count the atoms

    molecule.Frames().push_back(Frame(&molecule));
    Frame& frame=molecule.Frames().back();
    
    frame.XYZ().resize(molecule.Atoms().size());
    
    std::vector<Coordinate>::iterator ct=frame.XYZ().begin();
    while(std::getline(*m_file,line))
    {

      StringTokenizer token(line);
      size_t i;
      for(i=0;i<token.size();i++)
      {
        if(token[i]=="HETATM" || token[i]=="ATOM") //Do the stuff
        {
          int c=0;

          while(token[i+3+c].find(".")==std::string::npos)
          {
            c++;
          }

          ct->x()=std::atof(token.at(i+3+c).c_str());
          ct->y()=std::atof(token.at(i+4+c).c_str());
          ct->z()=std::atof(token.at(i+5+c).c_str());
          ++ct;
        }

      }
     if( line.find("MODEL") != std::string::npos /*|| line.find("END") != std::string::npos*/ ) break;
    }


  }

}

std::string PDBParser::ExtractAtomName(const std::string& s)
{
  PeriodicTable::const_iterator it;
  std::vector<PeriodicTable::const_iterator> atomnames;
  std::vector<size_t> cpos;
  for(it=GetPeriodicTable().begin();it!=GetPeriodicTable().end();++it)
  {
    size_t p=s.find(const_cast<std::string&>(it->first),0);
    if ( p == 0 )
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
