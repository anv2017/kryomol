/*****************************************************************************************
                            renderorbitals.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef RENDERORBITALS_H
#define RENDERORBITALS_H

#include <fstream>
#include <vector>
#include "mathtools.h"
#include "orbitalarray.h"
#include "grid.h"
#include "density.h"
#include "orbital.h"
#include "orbitaldata.h"
#include "transitionchange.h"
#include "renderexport.h"

namespace kryomol
{
class Frame;
class Coordinate;

class KRYOMOLRENDER_API RenderOrbitals
{
public:
    enum Axis {X,Y,Z};
    enum Basis {S, PX, PY, PZ, DXX, DXY, DXZ, DYY, DYZ, DZZ, DY0, DY1, DY2, DY3, DY4, FXXX, FXXY, FXXZ, FXYY, FXYZ, FXZZ, FYYY, FYYZ, FYZZ, FZZZ, FY0, FY1, FY2, FY3, FY4, FY5, FY6};
    RenderOrbitals() {}
    RenderOrbitals(Frame& frame);
    ~RenderOrbitals();

    void CalculateTotalDensity();
    void CalculateSpinDensity();
    void CalculateHomo();
    void CalculateLumo();
    void CalculateGridCoordinates();
    void CalculateAtomicOrbital (Basis basis, Coordinate& c, const OrbitalArray& exponential, OrbitalArray& atomicorbital);//(size_t ao, OrbitalArray& atomicorbital);
    void CalculateMolecularOrbital(size_t mo, OrbitalArray& molecularorbital, bool beta);
    void CalculateTransitionChange(size_t tc, OrbitalArray&  transition);
    void CalculateDensityChange(size_t tc, OrbitalArray&  transition);
    float CalculateContractionConstant(std::vector<float>& Xs, std::vector<float>& Alpha, Orbital::OrbitalType type);
    float CalculateAngularNormalizationConstant(Basis basis);
    void AdditionAtomicOrbital(const Coordinate& c, const OrbitalArray& atomicorbital, OrbitalArray& molecularorbital);

    void CleanMatrices(bool b_homo, bool b_lumo, bool b_total, bool b_spin, int b_orbital, int b_betaorbital);
    void ShowTotalDensity();
    void ShowSpinDensity();
    void ShowHomo();
    void ShowLumo();
    void ShowMolecularOrbital(size_t mo);
    void ShowBetaMolecularOrbital(size_t mo);
    void ShowTransitionChange(size_t tc);
    void ShowDensityChange(size_t dc);

    void SetThresholdAOSelector(float threshold) {m_thresholdAO = threshold; }
    void SetGridResolution(float resolution);
    void SetBeta(bool b) {m_beta = b;}

    Density& DensityData() { return m_density; }
    OrbitalData& OrbitalsData() { return m_orbitaldata; }
    float ThresholdAOSelector() { return m_thresholdAO;}
    float GridResolution() {return m_gridresolution;}
    Grid& GridData() {return m_grid;}

private:
    bool m_beta;
    bool m_diffuse;
    Coordinate m_centroid;
    Grid m_grid;
    Density m_density;
    OrbitalData m_orbitaldata;
    std::vector< std::vector<TransitionChange> > m_transitiondata;
    OrbitalArray  m_homo;
    OrbitalArray  m_lumo;
    OrbitalArray  m_totaldensity;
    OrbitalArray  m_spindensity;
    OrbitalArray m_xcoordinates;
    OrbitalArray m_ycoordinates;
    OrbitalArray m_zcoordinates;
    float m_thresholdAO;
    float m_gridresolution;
    std::vector< OrbitalArray > m_atomicorbitals;
    Fifo< int, OrbitalArray > m_molecularorbitals;
    Fifo< int, OrbitalArray > m_transitions;

};

}

#endif // RENDERORBITALS_H
