/*****************************************************************************************
                            molecule.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef MOLECULE_H
#define MOLECULE_H

#include <string>
#include <vector>
#include <list>



#include "atom.h"
#include "pdbtools.h"
#include "energy.h"
#include "mathtools.h"
#include "coreexport.h"
#include "frame.h"
#include "threshold.h"
#include "frequency.h"
/** */
namespace kryomol
{
  typedef std::vector<size_t> Ring;
  /** @brief Representation of a molecule

  In NMRDev molecules are multiconformational objects where each conformer is stored as a @see Molecule::Frame object
  */

  class MoleculePrivate;
  class FramePrivate;

  class KRYOMOLCORE_API Molecule
  {
    public:
      //enum Charge { ESP, MULLIKEN, NBO, AIM };
      struct pair
      {
      public:
        pair() {i=j=-1; value="";}
        pair(size_t a,size_t b,const std::string& s) {i=a;j=b;value=s;}
        size_t i; size_t j; std::string value;
      };
      struct triad
      {
      public:
        triad() { i=j=k=-1; value="";}
        triad(size_t a, size_t b, size_t c, const std::string& s) {i=a;j=b;k=c;value=s;}
        size_t i; size_t j; size_t k; std::string value;
      };
      struct quad
      {
      public:
        quad() { i=j=k=l=-1; value="";}
        quad(size_t a, size_t b, size_t c, size_t d,const std::string& s) {i=a;j=b;k=c;l=d;value=s;}
        size_t i; size_t j; size_t k; size_t l; std::string value;
      };

      class Coupling
      {
       public:
       Coupling(size_t i, size_t j, double value)  { _i=i; _j=j; _value=value; }
       size_t I() const { return _i; }
       size_t J() const { return _j; }
       double Value() const { return _value; }
       private:
       size_t _i;
       size_t _j;
       double _value;
      };
      /** Build an empty molecule*/
      Molecule();
      /** copy constructor*/
      Molecule ( const Molecule& mol );
      /** assignment operator*/
      Molecule& operator= ( const Molecule& mol );
      /** */
      ~Molecule();
      /** return frame i as a complete new molecule*/
      Molecule At ( size_t frame ) const;
      /** @return a const stl vector with all conformers*/
      const std::vector<Frame>& Frames() const;
      /** @return a stl vector with all conformers*/
      std::vector<Frame>& Frames();
      /** @return a const reference to the current conformer*/
      const Frame& CurrentFrame() const;
      /** @return a reference to the current conformer*/
      Frame& CurrentFrame();
      /** return the index of the current conformer*/
      size_t CurrentFrameIndex() const;
      /** set the active conformer to i*/
      void SetCurrentFrame ( size_t i ) ;
      /** @return a const stl vector of atoms in this molecule*/
      const std::vector<Atom>& Atoms() const;
      /** @return a const stl vector of atoms in this molecule*/
      std::vector<Atom>& Atoms() ;
      /** @return a const stl vector of bonds in this molecule*/
      const std::vector<Bond>& Bonds() const;
      /** @return a const stl vector of bonds in this molecule*/
      std::vector<Bond>& Bonds();
      /** @return the number of bonds between atoms i and j*/
      size_t NBonds(size_t i, size_t j) const;
      /** @return  populations*/
      const std::vector<double>& Populations() const;
      /** @return  populations*/
      std::vector<double>& Populations();
      /** returns molecular weight */
      double Weight() const { return CalculateWeight(); }
      /** move molecule to the centroid*/
      void MoveToCentroid();
      /** Get the rings of size size*/
      std::vector < std::vector<size_t> > Rings ( size_t size=0 ) const;
      /** obsolete */
      std::vector<int>  FindAtomByName ( const std::string& name ) const;
      /** Calculate Molecular Weight */
      double CalculateWeight() const;
      /** Get Hybridization type for atom i */
      Atom::hybridization Hybridization ( size_t atom ) const;
      /** automatically setup bonds based on atomic distances*/
      void SetBonds ( bool eachframe=false );
      /** return a list of PDB residues*/
      const std::vector<PDBResidue*>& Residues() const;
      /** return a list of PDB residues*/
      std::vector<PDBResidue*>& Residues();
      /** return a ordered list with the different element atomic numbers  present in the molecule*/
      std::list<int> Elements() const;
      /** return a ordered list with the different element atomic numbers  present in the molecule*/
      std::list<std::string> ElementSymbols() const;
      /** @return a vector with the indexes of the atoms connected to atom i*/
      std::vector<size_t> Neighbours ( size_t i ) const;
      /** return the index of the atom from its pdb name*/
      size_t IndexFromPDB(const std::string& pdbname,const std::string& resname,const std::string& resindex) const;
      /** Super impose frames to referance frame*/
      void SuperImpose(size_t refframe);
      /** superimpose all frames to reference frame f, minimizing the displacements for atoms*/
      void SuperImpose(size_t reframe, const std::vector<size_t>& atoms); 
      /** superimpose all frames to reference frame f, minimizing the displacements for atoms,using mass weighted coordinates*/
      void EckartTransform(size_t reframe, const std::vector<size_t>& atoms);
      std::vector<Coordinate>& InputOrientation() { return m_inputorientation; }
      const std::vector<Coordinate>& InputOrientation() const { return m_inputorientation; }
      std::vector< std::vector<Coordinate> >& GetMode() { return m_mode; }
      const std::vector< std::vector<Coordinate> >& GetMode() const { return m_mode; }
      //void SetMode(size_t index, const Coordinate c) { m_mode[index] = c; }

      void CalculateMassCenter(bool real=false);
      //void SetDihedral(size_t i, size_t j, size_t k, size_t l,float dihedral);
      void RotateBond(size_t i, size_t j, int sense);
      float GetAnglePlanes(std::vector<size_t>& v, bool degrees);
      std::vector < Coupling >& GetCouplings() { return m_couplings; }
      const std::vector < Coupling >& GetCouplings() const { return m_couplings; }
      void ResetSelection();
      unsigned int CountHeavyAtoms() const;

      std::string GetEnergyLevel() const;
      void SetEnergyLevel(std::string level) const;
      /** set a particular color for each conformation*/
      void SetColors();

  protected:
      void RotateNeighbours(size_t i, size_t j, int sense,double pass);
      Coordinate m_rotaxis[2];
      std::vector<Coordinate> m_inputorientation;
      std::vector < Coupling > m_couplings;
      std::vector<size_t> m_rotatedatoms;

      std::vector< std::vector<Coordinate> > m_mode;

  private:
      bool FindInShell(std::vector<size_t>& searched,size_t j) const;

  private:
      MoleculePrivate* m_private;

  };

  /** print molecule frames on stream s*/
  KRYOMOLCORE_API std::ostream& operator << ( std::ostream& s,const kryomol::Molecule& m );

}


#endif // !defined(AFX_MOLECULE_H__4A25AAA7_10C1_4F46_B82A_56A6C57A0AF3__INCLUDED_)
