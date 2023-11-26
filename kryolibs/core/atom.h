/*****************************************************************************************
                            atom.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef ATOM_H
#define ATOM_H


#include <string>
#include <vector>
#include <list>
#include <map>

#include "coordinate.h"
#include "coreexport.h"

/** @brief the global kryomol namespace
*/

namespace kryomol
{
  class PDBResidue;

  /** @brief representation of an atom

  This class stores basic info about an atom such its atomic number, mass, default color to be drawn and Pauling electronegativity
  */
  class KRYOMOLCORE_API Atom
  {
    public:
      class KRYOMOLCORE_API Isotope
      {
        /**
        @param mass, isotopic mass
        @param abundance, abundance (%)
        @param mult, multiplicity 2S+1
        @param lamba, magnetogyric ratio
        @param q, quadrupole
        */
        public:
        Isotope(int massnumber, double mass, double abundance, int mult, double lambda, double q) : _massnumber(massnumber), _mass(mass), _abundance(abundance), _mult(mult), _lambda(lambda), _q(q) {}
        int MassNumber() const { return _massnumber; }
        double Mass() const { return _mass; }
        double Abundance() const { return _abundance; }
        int Multiplicity() const { return _mult; }
        double MagnetoGyricRatio() const { return _lambda; }
        double Quadrupole() const { return _q; }
        private:
        /** mass number*/
        int _massnumber;
        /** isotopic mass in umas*/
        double _mass;
        /** abundance in %*/
        double _abundance;
        /** multiplicity (2S+1)*/
        int _mult;
        /** magnetogyric ratio*/
        double _lambda;
        /** quadrupole moment*/
        double _q;
      };
      
      /** @brief basic atomic data (Z, mass, electronegativity...)

      This class stores basic info about an atom such its atomic number, mass, default color to be drawn and Pauling electronegativity
      */
      class KRYOMOLCORE_API atomdata
      {
        public:

          /** RGB represenation of a color*/
          class rgb
          {
            public:
              rgb ( int lr, int lg, int lb ) : r ( lr ) , g ( lg ) ,b ( lb ) {}
              int r; int g; int b;
          };
          /** built atomdata

           This constructor is used in filling the periodic table
          @param z, atomic number
          @param m, atomic mass
          @param w, Van der Walss radium
          @param r, covalent atomic radium,
          @param el, Pauling Electronegativity
          @param c,  atom @see rgb color*/
          atomdata ( int z, double m, double w, double r, double el, rgb c,const std::list<Isotope>* idata=NULL ) : _Z ( z ), _mass ( m ), _vdw ( w ), _cr ( r ), _Pel ( el ), _color ( c ) , _idata(idata) {}
          /** atomic number*/
          int _Z;
          /** atomic mass*/
          double _mass;
          /** Van der Vals radium in Angstrom*/
          double _vdw;
          /** covalent atomic radium in Angstrom*/
          double _cr;
          /** Pauling electronegativity*/
          double _Pel;
          /** @see rgb color of the atom*/
          rgb _color;
          /** isotopic data*/
          const std::list<Isotope>* _idata;
      };

      /** hybridization of the atom*/
      enum hybridization { UNKNOWN=0, SP, SP2, SP3 };
      Atom();
      ~Atom();

      /** Construct atom with atomic symbol @param symbol
       
      This functions setup all atomic data adequately
      */
      Atom ( const std::string& symbol ) { SetSymbol ( symbol ); m_residue=NULL; }
      /** Construct atom with atomic number @param z

      This functions setup all atomic data adequately
      */
      Atom ( int z ) { SetZ ( z ); m_residue=NULL; }
      /** Return isotopically promediated atomic mass*/
      double AtomicMass() const { return m_data->second._mass; }
      /** fill @param r, @param g and @param b RGB parameters*/
      void Color ( float* r, float* g, float* b ) const;
      /** @return atomic number*/
      int Z() const { return m_data->second._Z; }
      /** @return Van der Waals radium*/
      double VdW() const { return m_data->second._vdw; }
      /** @return covalent radium*/
      double CovalentRadius() const { return m_data->second._cr; }
      /** @return Pauling Electronegativity*/
      double PaulingElectronegativity() const { return m_data->second._Pel; }
      /** @return the atomic symbol (uppercase)*/
      const std::string& Symbol() const { return m_data->first; }
      /** Set the atomic symbol to @param symbol and change Z */
      void SetSymbol ( const std::string& symbol );
      /** Set the atomic number to @param z and set symbol*/
      void SetZ ( int z );
      /** Set a pointer to the molecule pdb residues table*/
      void SetPDBResidue(PDBResidue* r) { m_residue=r; }
      /** return a pointer to the PDB residue*/
      PDBResidue* Residue() { return m_residue; }
      /** return a const pointer to the PDB residue*/
      const PDBResidue* Residue() const { return m_residue; }
      /** Set the PDB Name*/
      void SetPDBName(const std::string& name) { m_pdbname=name; }
      /** return the PDB name*/
      const std::string& PDBName() const { return m_pdbname; } 
      /** Get a reference for the active isotope*/
      const Isotope& CurrentIsotope() const;
      /**Set the active isotope to that of mass number m*/
      void SetCurrentIsotope(int m);
      /**return the list of atomic isotopes*/
      const std::list<Isotope>& Isotopes() const { return *m_data->second._idata; }
      bool IsSelected()   { return m_selected;}
      void Select(bool b) { m_selected=b;}

    protected:
      bool m_selected;

    private:
      double GetAtomicMass() const;
      const std::pair<const std::string,const atomdata>* m_data;
      const Isotope* m_isotope;
      PDBResidue* m_residue;
      std::string m_pdbname;

  };

  KRYOMOLCORE_API typedef std::map<const std::string,const kryomol::Atom::atomdata> PeriodicTable;

  /** @return a reference to the global periodic table*/
  KRYOMOLCORE_API const PeriodicTable& GetPeriodicTable();

  KRYOMOLCORE_API void BuildPeriodicTable();
  
  /** Return the atomic number Z for symbol */
  KRYOMOLCORE_API int ZFromSymbol(const std::string& symbol);
  
  /** Return the Symbol for element of atomic number z*/
   KRYOMOLCORE_API std::string SymbolFromZ(int z);
}

#endif // KRYOMOLATOM_H
