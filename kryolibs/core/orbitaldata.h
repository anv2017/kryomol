/*****************************************************************************************
                            orbitaldata.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef ORBITALDATA_H
#define ORBITALDATA_H

#include <iostream>
#include <vector>
#include "mathtools.h"
#include "coreexport.h"
#include "orbital.h"
#include "basiscenter.h"

namespace kryomol
{

/** @brief atomic and molecular orbitals data

This class manage atomic and molecular orbitals data from compiled from parser file*/
class KRYOMOLCORE_API OrbitalData
{
public:
    OrbitalData();
    ~OrbitalData() {}

    int Homo() { return m_homo; }
    int Lumo() { return m_lumo; }
    int TypeD() { return m_typeD; }
    int TypeF() { return m_typeF; }
    const D2Array<float>& Coefficients() const { return m_coefficients;}
    D2Array<float>& Coefficients() { return m_coefficients;}
    const std::vector<float>& Eigenvalues() const { return m_eigenvalues; }
    std::vector<float>& Eigenvalues() { return m_eigenvalues; }
    const std::vector<float>& Occupations() const { return m_occupations; }
    std::vector<float>& Occupations() { return m_occupations; }
    const D2Array<float>& BetaCoefficients() const { return m_betacoefficients;}
    D2Array<float>& BetaCoefficients() { return m_betacoefficients;}
    const std::vector<float>& BetaEigenvalues() const { return m_eigenvalues; }
    std::vector<float>& BetaEigenvalues() { return m_betaeigenvalues; }
    const std::vector<Orbital>& Orbitals() const { return m_orbitals; }
    std::vector<Orbital>& Orbitals() { return m_orbitals; }
    const std::vector<BasisCenter>& BasisCenters() const { return m_basiscenters; }
    std::vector<BasisCenter>& BasisCenters() { return m_basiscenters; }

    void SetHomo(int homo) { m_homo = homo; }
    void SetLumo(int lumo) { m_lumo = lumo; }
    void SetTypeD(int typeD) { m_typeD = typeD; }
    void SetTypeF(int typeF) { m_typeF = typeF; }
    void SetCoefficients(D2Array<float> matrix) { m_coefficients = matrix; }
    void SetEigenvalues(std::vector<float> eigenvalues) { m_eigenvalues = eigenvalues; }
    void SetOrbitals(std::vector<Orbital> orbitals) { m_orbitals = orbitals; }    
    void SetBetaCoefficients(D2Array<float> matrix) { m_betacoefficients = matrix; }
    void SetBetaEigenvalues(std::vector<float> eigenvalues) { m_betaeigenvalues = eigenvalues; }
    void SetBasisCenters(std::vector<BasisCenter> basis) {m_basiscenters = basis;}
    void SetOccupations(const std::vector<float>& v) { m_occupations=v; }

private:
    int m_homo;
    int m_lumo;
    int m_typeD;
    int m_typeF;
    D2Array<float> m_coefficients;
    std::vector<float> m_eigenvalues;
    std::vector<float> m_occupations;
    D2Array<float> m_betacoefficients;
    std::vector<float> m_betaeigenvalues;
    std::vector<Orbital> m_orbitals;
    std::vector<BasisCenter> m_basiscenters;

};

}

#endif // ORBITALDATA_H
