/*****************************************************************************************
                            atom.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "atom.h"
#include "stringtools.h"

using namespace kryomol;

std::map<const std::string,const Atom::atomdata> m_periodictable;
std::vector< std::list< Atom::Isotope> > m_nmrtable;

/** Buld an undefined atom*/
Atom::Atom() : m_data ( NULL ), m_isotope(NULL), m_residue(NULL)
{
    m_selected=false;
}

Atom::~Atom()
{
}

/** Set atomica number for the atom, symbol and chemical properties
will be updated*/
void Atom::SetZ ( int z )
{
  std::map<const std::string,const Atom::atomdata>::const_iterator it;
  for ( it=GetPeriodicTable().begin();it!=GetPeriodicTable().end();++it )
  {
    if ( ( *it ).second._Z == z )
    {
      m_data=& ( *it );
      m_isotope=NULL;
      return;
    }
  }
}

/** Set current isotope to that of mass number m*/
void Atom::SetCurrentIsotope(int m)
{
  std::list<Isotope>::const_iterator it;
  for(it=m_data->second._idata->begin();it!=m_data->second._idata->end();++it)
  {
    if (it->MassNumber() == m )
      m_isotope=&(*it);
  }
}

/** return the current isotope for the atom*/
const Atom::Isotope& Atom::CurrentIsotope() const
{
  if ( m_isotope == NULL )
     return m_data->second._idata->front();
  else
  return *m_isotope;
} 

/** Set chemical symbol for the atom and associated chemical data*/
void Atom::SetSymbol ( const std::string& symbol )
{
  std::map<const std::string,const Atom::atomdata>::const_iterator it;
  it = m_periodictable.find(symbol);
  if ( it == m_periodictable.end() )
    it=m_periodictable.find("X");
  
  m_data=&(*it);
  m_isotope=NULL;
}

/** Convert 0-255 RGB values to 0-1 float rgb values*/
void Atom::Color ( float* r, float* g, float* b ) const
{
  *r=m_data->second._color.r/255.;
  *g=m_data->second._color.g/255.;
  *b=m_data->second._color.b/255.;
}

const float nodata=-1;

/** \brief buld the nmr periodic table

Bult a stl map based NMR periodic table
Each entry is a std::pair< std::string, data> element where the first 
element is the chemical symbol of the element and the second one is a
atomdata structure. The atomdata structure stores a pointer to the NMRTable
with the isotopic information.
This function should be called at the beginning of any qryomol based application
Dummy atoms are represented by the symbol X and have 0 as atomic and mass numbers
*/
void kryomol::BuildPeriodicTable()
{
  typedef Atom::atomdata data;
  typedef Atom::Isotope idata;
  
  std::list<Atom::Isotope> tdata;
  
  //hydrogen
  tdata.push_back( idata(1,1.007825032,99.9885,2,26.7522128,0));
  tdata.push_back( idata(2,2.014101778,0.0115,3,4.10662791,2.8e-3));
  tdata.push_back( idata(3,3.0160492675,0,2,28.5349779,0));
  m_nmrtable.push_back(tdata);
  m_periodictable.insert ( std::pair<std::string,data> ( "H",data ( 1,1.00794,1.20,0.370,2.20,data::rgb ( 170,170,170 ), &m_nmrtable.back() ) ) );
  
  
  //helium
  m_periodictable.insert ( std::pair<std::string,data> ( "He",data ( 2,4.002602,1.40,0.320,nodata,data::rgb ( 233,194,243 ) ) ) );
  
  //lithium
  tdata.clear();
  tdata.push_back( idata(6,6.015122,7.59,3,3.9371709,-0.0808) );
  tdata.push_back( idata(7,7.016004,92.41,4,10.3977013,-4.01) );
  m_nmrtable.push_back(tdata);
  m_periodictable.insert ( std::pair<std::string,data> ( "Li",data ( 3,6.941,1.82,1.34,0.98,data::rgb ( 50,73,255 ),&m_nmrtable.back() ) ) );
  
  //berillium
  m_periodictable.insert ( std::pair<std::string,data> ( "Be",data ( 4,9.012182,2*0.90,0.90,1.57,data::rgb ( 169,170,255 ) ) ) );
  //boron
  m_periodictable.insert ( std::pair<std::string,data> ( "B",data ( 5,10.811,2*0.82,0.82,2.04,data::rgb ( 193,211,201 ) ) ) );
  
  //carbon
  tdata.clear();
  tdata.push_back( idata(13,13.003354,1.07,2,6.728284,0)) ;
  m_nmrtable.push_back(tdata);
  m_periodictable.insert ( std::pair<std::string,data> ( "C",data ( 6,12.0107,1.70,0.77,2.55,data::rgb ( 050,200,050 ),&m_nmrtable.back() ) ) );
  
  //nitrogen
  tdata.clear();
  tdata.push_back( idata(15,15.000108,0.368,2,-2.71261804,0) );
  tdata.push_back( idata(14,14.003074,99.632,3,1.9337792,0) );
  m_nmrtable.push_back(tdata);
  m_periodictable.insert ( std::pair<std::string,data> ( "N",data ( 7,14.0067,1.55,0.75,3.04,data::rgb ( 050,050,200 ),&m_nmrtable.back() ) ) );
  
  //oxygen
  m_periodictable.insert ( std::pair<std::string,data> ( "O",data ( 8,15.9994,1.52,0.73,3.44,data::rgb ( 200,050,050 ) ) ) );
  
  //fluorine
  tdata.clear();
  tdata.push_back( idata(19,18.99840322,100.,2,25.18148,0) );
  m_nmrtable.push_back(tdata);
  m_periodictable.insert ( std::pair<std::string,data> ( "F",data ( 9,18.9984032,1.47,0.71,3.98,data::rgb ( 181,243,184 ),&m_nmrtable.back() ) ) );
  
  m_periodictable.insert ( std::pair<std::string,data> ( "Ne",data ( 10,20.1797,1.54,0.69,nodata,data::rgb ( 233,194,243 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Na",data ( 11,22.989770,2.27,1.54,0.93,data::rgb ( 50,73,255 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Mg",data ( 12,24.3050,0,0,1.31,data::rgb ( 169,170,255 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Al",data ( 13,26.981538,2*1.18,1.18,1.61,data::rgb ( 193,211,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Si",data ( 14,28.0855,0,0,1.90,data::rgb ( 193,211,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "P",data ( 15,30.973761,1.80,1.06,2.19,data::rgb ( 197,117,243 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "S",data ( 16,32.065,1.80,1.02,2.58,data::rgb ( 228,243,91 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Cl",data ( 17,35.453,0,0,3.16,data::rgb ( 181,243,184 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Ar",data ( 18,39.948 ,0,0,nodata,data::rgb ( 128,128,128 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "K",data ( 19,39.0983,0,0,0.82,data::rgb ( 50,73,255 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Ca",data ( 20,40.078,0,0,1.00,data::rgb ( 169,170,255 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Sc",data ( 21,44.955910,0,0,1.36,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Ti",data ( 22,47.867,0,0,1.54,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "V",data ( 23,50.9415,0,0,1.63,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Cr",data ( 24,51.9961,0,0,1.66,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Mn",data ( 25,54.938049,0,0,1.55,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Fe",data ( 26,55.845,0,0,1.03,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Co",data ( 27,58.933200,0,0,1.88,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Ni",data ( 28,58.6934,0,0,1.91,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Cu",data ( 29,63.546,0,0,1.90,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Zn",data ( 30,65.409,0,0,1.65,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Ga",data ( 31,69.723,0,0,1.81,data::rgb ( 193,211,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Ge",data ( 32,72.64,0,0,2.01,data::rgb ( 193,211,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "As",data ( 33,74.92160,0,0,2.18,data::rgb ( 197,117,243 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Se",data ( 34,78.96,0,0,2.55,data::rgb ( 228,243,91 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Br",data ( 35,79.904,0,0,2.96,data::rgb ( 181,243,184 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Kr",data ( 36,83.798,0,0,3.00,data::rgb ( 233,194,243 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Rb",data ( 37,85.4678,0,0,0.82,data::rgb ( 50,73,255 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Sr",data ( 38,87.62,0,0,0.95,data::rgb ( 169,170,255 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Y",data ( 39,88.90585,0,0,1.22,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Zr",data ( 40,91.224,0,0,1.33,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Nb",data ( 41,92.90638,0,0,1.6,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Mo",data ( 42,95.94,0,0,2.16,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Tc",data ( 43,98.,0,0,1.9,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Ru",data ( 44,101.07,0,0,2.2,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Rh",data ( 45,102.90550,0,0,2.28,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Pd",data ( 46,106.42,0,0,2.20,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Ag",data ( 47,107.8682,0,0,1.93,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Cd",data ( 48,112.411,0,0,1.69,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "In",data ( 49,114.818,0,0,1.78,data::rgb ( 193,211,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Sn",data ( 50,118.710,0,0,1.96,data::rgb ( 193,211,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Sb",data ( 51,121.760,0,0,2.05,data::rgb ( 197,117,243 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Te",data ( 52,127.60,0,0,2.1,data::rgb ( 228,243,91 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "I",data ( 53,126.90447,0,0,2.66,data::rgb ( 181,243,184 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Xe",data ( 54,131.293,0,0,2.6,data::rgb ( 233,194,243 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Cs",data ( 55,132.90545,0,0,0.79,data::rgb ( 50,73,255 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Ba",data ( 56,137.327,0,0,0.89,data::rgb ( 169,170,255 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "La",data ( 57,138.9055,0,0,1.10,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Ce",data ( 58,140.116,0,0,1.12,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Pr",data ( 59,140.90765,0,0,1.13,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Nd",data ( 60,144.24,0,0,1.14,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Pm",data ( 61,145. ,0,0,nodata,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Sm",data ( 62,150.36,0,0,1.17,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Eu",data ( 63,151.964 ,0,0,nodata,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Gd",data ( 64,157.25,0,0,1.20,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Tb",data ( 65,158.92534 ,0,0,nodata,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Dy",data ( 66,162.500,0,0,1.22,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Ho",data ( 67,164.93032,0,0,1.13,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Er",data ( 68,167.259,0,0,1.24,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Tm",data ( 69,168.932421,0,0,1.25,data::rgb ( 100,100,100 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Yb",data ( 70,173.04 ,0,0,nodata,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Lu",data ( 71,174.967,0,0,1.27,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Hf",data ( 72,178.49,0,0,1.3,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Ta",data ( 73,180.9479,0,0,1.5,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "W",data ( 74,183.84,0,0,2.36,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Re",data ( 75,186.207,0,0,1.9,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Os",data ( 76,190.23,0,0,2.2,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Ir",data ( 77,192.217,0,0,2.20,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Pt",data ( 78,195.078,0,0,2.28,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Au",data ( 79,196.96655,0,0,2.54,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Hg",data ( 80,200.59,0,0,2.00,data::rgb ( 214,41,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Tl",data ( 81,204.3833,0,0,1.62,data::rgb ( 193,211,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Pb",data ( 82,207.2,0,0,2.33,data::rgb ( 193,211,201 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Bi",data ( 83,208.98038,0,0,2.02,data::rgb ( 197,117,243 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Po",data ( 84,209.,0,0,2.0,data::rgb ( 228,243,91 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "At",data ( 85,210.,0,0,2.2,data::rgb ( 181,243,184 ) ) ) );
  m_periodictable.insert ( std::pair<std::string,data> ( "Rn",data ( 86,222. ,0,0,nodata,data::rgb ( 233,194,243 ) ) ) );
  //dummy atom
  m_periodictable.insert ( std::pair<std::string,data> ( "X",data (0,0. ,0,0,nodata,data::rgb ( 255,0,0 ) ) ) );



}

const kryomol::PeriodicTable& kryomol::GetPeriodicTable()
{
  return m_periodictable;
}

int kryomol::ZFromSymbol(const std::string& symbol)
{
  int z=0;
  const kryomol::PeriodicTable& table=GetPeriodicTable();
  std::map<const std::string,const kryomol::Atom::atomdata>::const_iterator it;
  for ( it=table.begin();it!=table.end();++it )
  {
    if ( ( *it ).first == symbol )
    {
      z = (*it).second._Z;
      return z;
    }
  }
  return z;
}

std::string kryomol::SymbolFromZ(int z)
{
  std::string s;
  const PeriodicTable& table=GetPeriodicTable();
  std::map<const std::string,const Atom::atomdata>::const_iterator it;
  for ( it=table.begin();it!=table.end();++it )
  {
    if ( ( *it ).second._Z == z)
    {
        s = (*it).first;
        return s;
    }
  }
  return s;
} 
