/*****************************************************************************************
                            frame.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef FRAME_H
#define FRAME_H

#include <vector>
#include "coordinate.h"
#include "atom.h"
#include "bond.h"
#include "energy.h"
#include "coreexport.h"
#include "frequency.h"
#include "threshold.h"
#include "mathtools.h"
#include "orbitaldata.h"
#include "transitionchange.h"
#include "electronicdensity.h"
#include "renderdensity.h"
#include "grid.h"

namespace kryomol
{
  class Molecule;
  class FramePrivate;

  /** @brief Representation of a conformer*/
  class KRYOMOLCORE_API Frame
  {
    public:
      enum Charge { ESP, MULLIKEN, NBO, AIM };
      Frame ( Molecule* molecule );
      Frame ( const Frame& frame );
      Frame& operator = ( const Frame& frame );
      ~Frame();
      /** @return a const vector of bonds*/
      const std::vector<Bond>& Bonds() const;
      /** @return a const vector of bonds*/
      std::vector<Bond>& Bonds();
      /** @return coordinates for this frame*/
      const std::vector<Coordinate>& XYZ() const ;
      /** @return coordinates for this frame*/
      std::vector<Coordinate>& XYZ();
      /** @return the inertia tensor*/
      D2Array<double> InertiaTensor() const;
      /** @return the Gyration tensor*/
      D2Array<double> GyrationTensor() const;
      /** set the dihedral angle in radians */
      void SetDihedral ( size_t i,size_t j, size_t k, size_t l, float dihedral );
      /** rotate bond between atoms i and j*/
      void RotateBond ( size_t i,size_t j,float angle, bool clockwise );
      /** @return the kinetic energy of this frame*/
      Energy& KineticEnergy() ;
      /** @return the kinetic energy of this frame*/
      const Energy& KineticEnergy() const;
      /** @return the potential energy of this frame*/
      const Energy& PotentialEnergy() const;
      /** @return the potential energy of this frame*/
      Energy& PotentialEnergy() ;
      /** @return total energy as the sum of potential and kinetic energy*/
      Energy TotalEnergy() const ;
      /** @return the RMS force of this frame*/
      const Energy& RMSForce() const;
      /** @return the RMS force of this frame*/
      Energy& RMSForce() ;
      /** @return the distance in Angstroms between atoms i and j*/
      float Distance ( int i, int j ) const;
      /** @return angle between atoms i j and l*/
      float Angle ( int i, int j, int l ) const;
      /** @return dihedral angle between atoms i,j,k,l*/
      float Dihedral ( int i, int j, int k, int l ) const;
      /** @return position of the Center of mass*/
      Coordinate MassCenter() const;
      /** @return position of the centroid*/
      Coordinate Centroid() const;
      /** @return rings of size @param size or all rings if 0*/
      std::vector < std::vector<size_t> > Rings ( size_t size=0 ) const;
      /** @return 3D vector betwwen atoms i and j*/
      Coordinate Vector ( int i, int j ) const;
      /** @return a const pointer to the parent molecule of this conformer*/
      const Molecule* ParentMolecule() const { return m_molecule; }
      /** @return a pointer to the parent molecule of this conformer*/
      Molecule* ParentMolecule() { return m_molecule; }
      /** @return a vector with the chemical shift tensors*/
      const std::vector< D2Array<double> >& CShiftTensors() const;
      /** @return a vector with the chemical shift tensors*/
      std::vector< D2Array<double> >& CShiftTensors();
      /** @return the gradient coordinates for this frame*/
      std::vector<Coordinate>& Gradient();
      /** @return the gradient coordinates for this frame*/
      const std::vector<Coordinate>& Gradient() const;
      /** @return the orbital data of this frame*/
      OrbitalData& OrbitalsData() ;
      /** @return the orbitaldata of this frame*/
      const OrbitalData& OrbitalsData() const;
      /** @return the transitionchanges of this frame*/
      std::vector< std::vector<TransitionChange> >& TransitionChanges() ;
      /** @return the transitionchanges of this frame*/
      const std::vector< std::vector<TransitionChange> >& TransitionChanges() const;
      /** @return the electronic density data of this frame*/
      ElectronicDensity& ElectronicDensityData() ;
      /** @return the electronic density data of this frame*/
      const ElectronicDensity& ElectronicDensityData() const;
      /** @return the vector with the positive electronic density data of this frame*/
      std::vector<RenderDensity>& PositiveDensity() ;
      /** @return the vector with the positive electronic density data of this frame*/
      const std::vector<RenderDensity>& PositiveDensity() const;
      /** @return the vector with the negative electronic density data of this frame*/
      std::vector<RenderDensity>& NegativeDensity() ;
      /** @return the vector with the negative electronic density data of this frame*/
      const std::vector<RenderDensity>& NegativeDensity() const;
      /** @return a grid for calculating the molecular orbitals*/
      Grid& GetGrid();
      /** @return a grid for calculating the molecular orbitals*/
      const Grid& GetGrid() const;
      /** @return the coordinates (xmin,ymin,zmin) (xmax,ymax,zmax) for the box containing the molecule*/
      std::pair<Coordinate,Coordinate> Box() const;
      Coordinate Dipole();
      D1Array<double>& GetForces() ;
      D2Array<double>& GetHessian()  ;
      /** get the HSL representative color for the conformer*/
      void GetColor(float& h,float& s,float& l) const;
      /** set the HSL representative color for the conformer*/
      void SetColor(float h,float s,float l);
      void SetEnergy(double e, const std::string& level="");
      void SetS2(double e) ;
      void SetRMSForce(double f) ;
      void SetMaximumForce(double f) ;
      void SetRMSDisplacement(double d) ;
      void SetMaximumDisplacement(double d);
      void SetDipole(const Coordinate& c);
      void SetGrid(const Grid& grid);
      void SetOrbitalData(const OrbitalData& orbitaldata);
      void SetTransitionChanges(const std::vector< std::vector<TransitionChange> >& transitions);
      void SetElectronicDensityData(const ElectronicDensity& electronicdensity);
      void SetPositiveDensity(const std::vector<RenderDensity>& positivedensity);
      void SetNegativeDensity(const std::vector<RenderDensity>& negativedensity);

      double GetEnergy() const;
      double GetRMSForce() const ;
      double GetMaximumForce() const ;
      double GetRMSDisplacement() const ;
      double GetMaximumDisplacement() const ;
      double GetS2() const ;
      Threshold GetThreshold() const ;
      void SetThreshold(const Threshold& thr) ;
      const std::vector<double>* Charges( Charge type) const;
      std::vector<double>* Charges( Charge type);
      std::vector<Frequency> &GetFrequencies() ;
      const std::vector<Frequency>& GetFrequencies() const;
      std::vector<Spectralline>& GetSpectralLines();
      const std::vector<Spectralline>& GetSpectralLines() const;
      Coordinate GetDipole() const;

      void AllocateVectors(bool rottransmodes=false);
      void AllocateForces();
      void AllocateHessian();

      void CalculateGrid(float resolution);

    private:
      void RotateNeighbours ( const Coordinate& axisorigin,const Coordinate& axisend,size_t i,size_t j,float pass,bool clockwise,std::vector<size_t>& rotatedatoms );

    private:
      Molecule* m_molecule;
      FramePrivate* m_private;

  };

  /** print frame coordinates on stream s*/
  KRYOMOLCORE_API std::ostream& operator << ( std::ostream& s,const kryomol::Frame& m );
}

#endif
