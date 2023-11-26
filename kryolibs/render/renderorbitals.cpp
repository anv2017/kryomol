/*****************************************************************************************
                            renderorbitals.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <QApplication>
#include <QCursor>

#ifdef WITH_TIMERS
#include <QTime>
#endif

#include <QDebug>


#include "renderorbitals.h"
#include "frame.h"
#include "grid.h"
#include "stringtools.h"
#include "physicalconstants.h"
#include "orbitalarray.h"


using namespace kryomol;

RenderOrbitals::RenderOrbitals(Frame& frame)
{
    if (!frame.OrbitalsData().BasisCenters().empty())
    {
        m_thresholdAO = 0.001;
        m_gridresolution = 0.2;

        m_orbitaldata = frame.OrbitalsData();
        m_transitiondata = frame.TransitionChanges();

        m_centroid = frame.Centroid();
        qDebug() << "CENTROIDE: " << m_centroid.x() << m_centroid.y() << m_centroid.z() << "******************"  << endl;

        frame.CalculateGrid(m_gridresolution);
        m_grid = frame.GetGrid();
        m_density = Density(m_grid.Nx(), m_grid.Ny(), m_grid.Nz(), m_grid.Nl(), m_grid.Step(), m_grid.Step(), m_grid.Step(), m_grid.Step(), Coordinate(-m_grid.X()/2,-m_grid.Y()/2,-m_grid.Z()/2));
    }
    else
    {
        m_density = Density(frame.ElectronicDensityData().Nx(), frame.ElectronicDensityData().Ny(), frame.ElectronicDensityData().Nz(), frame.ElectronicDensityData().Dx(), frame.ElectronicDensityData().Dy(), frame.ElectronicDensityData().Dz(), frame.ElectronicDensityData().Origin());
        m_totaldensity = frame.ElectronicDensityData().Density();
    }
}

RenderOrbitals::~RenderOrbitals()
{}


void RenderOrbitals::SetGridResolution(float resolution)
{
    m_gridresolution = resolution;
    m_grid.SetStep(m_gridresolution);

    m_density = Density(m_grid.Nx(), m_grid.Ny(), m_grid.Nz(), m_grid.Nl(), resolution, resolution, resolution, resolution, Coordinate(-m_grid.X()/2,-m_grid.Y()/2,-m_grid.Z()/2));

}

void RenderOrbitals::CalculateGridCoordinates()
{
    size_t nx = m_density.Nx();
    size_t ny = m_density.Ny();
    size_t nz = m_density.Nz();

    m_xcoordinates = OrbitalArray(nx,ny,nz);
    m_ycoordinates = OrbitalArray(nx,ny,nz);
    m_zcoordinates = OrbitalArray(nx,ny,nz);

    float x= m_density.Origin().x();
    float y= m_density.Origin().y();
    float z= m_density.Origin().z();

    size_t i=0;
    while (i<nx)
    {
        size_t j=0;
        while (j<ny)
        {
            size_t k=0;
            while (k<nz)
            {
                m_xcoordinates(i,j,k) = x;
                m_ycoordinates(i,j,k) = y;
                m_zcoordinates(i,j,k) = z;

                ++k;
                z += m_grid.Step();
            }
            ++j;
            z= m_density.Origin().z();
            y += m_grid.Step();
        }
        ++i;
        y= m_density.Origin().y();
        x += m_grid.Step();
    }
}

void RenderOrbitals::CalculateAtomicOrbital(Basis basis,Coordinate &coordinate, const OrbitalArray &exponential, OrbitalArray &atomicorbital)
{
    switch (basis)
    {
        case RenderOrbitals::S:
            atomicorbital.CalculateOrbitalS(exponential);
            break;
        case RenderOrbitals::PX:
            atomicorbital.CalculateOrbitalPx(exponential,coordinate.x(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::PY:
            atomicorbital.CalculateOrbitalPy(exponential,coordinate.y(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::PZ:
            atomicorbital.CalculateOrbitalPz(exponential,coordinate.z(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::DXX:
            atomicorbital.CalculateOrbitalDxx(exponential,coordinate.x(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::DYY:
            atomicorbital.CalculateOrbitalDyy(exponential,coordinate.y(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::DZZ:
            atomicorbital.CalculateOrbitalDzz(exponential,coordinate.z(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::DXY:
            atomicorbital.CalculateOrbitalDxy(exponential,coordinate.x(),coordinate.y(),m_grid.L(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::DXZ:
            atomicorbital.CalculateOrbitalDxz(exponential,coordinate.x(),coordinate.z(),m_grid.L(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::DYZ:
            atomicorbital.CalculateOrbitalDyz(exponential,coordinate.y(),coordinate.z(),m_grid.L(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::FXXX:
            atomicorbital.CalculateOrbitalFxxx(exponential,coordinate.x(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::FYYY:
            atomicorbital.CalculateOrbitalFyyy(exponential,coordinate.y(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::FZZZ:
            atomicorbital.CalculateOrbitalFzzz(exponential,coordinate.z(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::FXXY:
            atomicorbital.CalculateOrbitalFxxy(exponential,coordinate.x(),coordinate.y(),m_grid.L(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::FXYY:
            atomicorbital.CalculateOrbitalFxyy(exponential,coordinate.x(),coordinate.y(),m_grid.L(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::FXXZ:
            atomicorbital.CalculateOrbitalFxxz(exponential,coordinate.x(),coordinate.z(),m_grid.L(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::FXZZ:
            atomicorbital.CalculateOrbitalFxzz(exponential,coordinate.x(),coordinate.z(),m_grid.L(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::FYYZ:
            atomicorbital.CalculateOrbitalFyyz(exponential,coordinate.y(),coordinate.z(),m_grid.L(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::FYZZ:
            atomicorbital.CalculateOrbitalFyzz(exponential,coordinate.y(),coordinate.z(),m_grid.L(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::FXYZ:
            atomicorbital.CalculateOrbitalFxyz(exponential,coordinate.x(),coordinate.y(),coordinate.z(),m_grid.L(),m_grid.L(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::DY0:
            atomicorbital.CalculateOrbitalDY0(exponential,coordinate.x(),coordinate.y(),coordinate.z(),m_grid.L(),m_grid.L(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::DY1:
            atomicorbital.CalculateOrbitalDY1(exponential,coordinate.x(),coordinate.z(),m_grid.L(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::DY2:
            atomicorbital.CalculateOrbitalDY2(exponential,coordinate.y(),coordinate.z(),m_grid.L(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::DY3:
            atomicorbital.CalculateOrbitalDY3(exponential,coordinate.x(),coordinate.y(),m_grid.L(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::DY4:
            atomicorbital.CalculateOrbitalDY4(exponential,coordinate.x(),coordinate.y(),m_grid.L(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::FY0:
            atomicorbital.CalculateOrbitalFY0(exponential,coordinate.x(),coordinate.y(),coordinate.z(),m_grid.L(),m_grid.L(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::FY1:
            atomicorbital.CalculateOrbitalFY1(exponential,coordinate.x(),coordinate.y(),coordinate.z(),m_grid.L(),m_grid.L(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::FY2:
            atomicorbital.CalculateOrbitalFY2(exponential,coordinate.x(),coordinate.y(),coordinate.z(),m_grid.L(),m_grid.L(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::FY3:
            atomicorbital.CalculateOrbitalFY3(exponential,coordinate.x(),coordinate.y(),coordinate.z(),m_grid.L(),m_grid.L(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::FY4:
            atomicorbital.CalculateOrbitalFY4(exponential,coordinate.x(),coordinate.y(),coordinate.z(),m_grid.L(),m_grid.L(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::FY5:
            atomicorbital.CalculateOrbitalFY5(exponential,coordinate.x(),coordinate.y(),m_grid.L(),m_grid.L(),m_grid.Step());
            break;
        case RenderOrbitals::FY6:
            atomicorbital.CalculateOrbitalFY6(exponential,coordinate.x(),coordinate.y(),m_grid.L(),m_grid.L(),m_grid.Step());
            break;

        default:
            break;


        }
}


void RenderOrbitals::CalculateMolecularOrbital(size_t mo, OrbitalArray& molecularorbital, bool beta)
{
#ifdef WITH_TIMERS
    QTime timer;
    timer.start();
#endif
    QApplication::setOverrideCursor(QCursor(Qt::BusyCursor));

    molecularorbital.SetToZero();

    D2Array<float> coefficients;
    if (beta)
        coefficients = m_orbitaldata.BetaCoefficients();
    else
        coefficients = m_orbitaldata.Coefficients();

    float coeff_j_i;
    size_t N = 0;
    for (size_t j=0; j<m_orbitaldata.BasisCenters().size(); ++j)
    {
        Coordinate atom = m_orbitaldata.BasisCenters().at(j).Atom();
        Coordinate centersubgrid = Coordinate(m_grid.Step()*floor((atom.x()-m_density.Origin().x())/m_grid.Step())+m_density.Origin().x(), m_grid.Step()*floor((atom.y()-m_density.Origin().y())/m_grid.Step())+m_density.Origin().y(), m_grid.Step()*floor((atom.z()-m_density.Origin().z())/m_grid.Step())+m_density.Origin().z());
        Coordinate c = atom-centersubgrid;
        /*qDebug() << "STEP: " << m_grid.Step() << endl;
        qDebug() << "ORIGIN: " << m_density.Origin().x() << m_density.Origin().y() << m_density.Origin().z() << endl;
        qDebug() << "GRID: " << m_grid.X() << m_grid.Y() << m_grid.Z() << endl;
        qDebug() << "ATOM: " << atom.x() << atom.y() << atom.z() << endl;
        qDebug() << "CENTER SUBGRID: " << centersubgrid.x() << centersubgrid.y() << centersubgrid.z() << endl;
        qDebug() << "SUBGRID CALCULATION:" << floor((atom.x()-m_density.Origin().x())/m_grid.Step()) << m_grid.Step()*floor((atom.x()-m_density.Origin().x())/m_grid.Step()) << m_density.Origin().x() << floor((atom.y()-m_density.Origin().y())/m_grid.Step()) << m_grid.Step()*floor((atom.y()-m_density.Origin().y())/m_grid.Step()) << m_density.Origin().y() << floor((atom.z()-m_density.Origin().z())/m_grid.Step()) << m_grid.Step()*floor((atom.z()-m_density.Origin().z())/m_grid.Step()) << m_density.Origin().z() << endl;
*/

        for (size_t i=0; i<m_orbitaldata.BasisCenters().at(j).Orbitals().size(); ++i)
        {            
            m_grid.CalculateSubgrid(m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha());            

            OrbitalArray atomicorbital(m_grid.Nl(),m_grid.Nl(),m_grid.Nl());
            OrbitalArray exponential(m_grid.Nl(),m_grid.Nl(),m_grid.Nl());
            OrbitalArray exponentialP(m_grid.Nl(),m_grid.Nl(),m_grid.Nl());

            if (m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Type() == Orbital::S)
            {
                coeff_j_i = coefficients(N,mo);

                if (fabs(coeff_j_i) > m_thresholdAO)
                {
                    atomicorbital.SetToZero();
                    for (size_t k=0; k<m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha().size(); ++k)
                    {
                        atomicorbital.CalculateExponential(m_grid.L(),m_grid.L(),m_grid.L(),m_grid.Step(),c.x(),c.y(),c.z(),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha().at(k),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Xs().at(k));
                    }
                    float A = CalculateContractionConstant(m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Xs(),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha(),Orbital::S);
                    coeff_j_i*=A;
                    atomicorbital*=coeff_j_i;
                    AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                }
                N+=1;
            }

            if (m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Type() == Orbital::SP)
            {
                if ( (fabs(coefficients(N,mo)) > m_thresholdAO) || (fabs(coefficients(N+1,mo)) > m_thresholdAO) || (fabs(coefficients(N+2,mo)) > m_thresholdAO) || (fabs(coefficients(N+3,mo)) > m_thresholdAO) )
                {
                    exponential.SetToZero();
                    exponentialP.SetToZero();
                    for (size_t k=0; k<m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha().size(); ++k)
                    {
                        exponential.CalculateExponential(m_grid.L(),m_grid.L(),m_grid.L(),m_grid.Step(),c.x(),c.y(),c.z(),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha().at(k),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Xs().at(k));
                        exponentialP.CalculateExponential(m_grid.L(),m_grid.L(),m_grid.L(),m_grid.Step(),c.x(),c.y(),c.z(),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha().at(k),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Xp().at(k));
                    }
                    float A = CalculateContractionConstant(m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Xs(),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha(),Orbital::S);
                    float Ap = CalculateContractionConstant(m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Xp(),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha(),Orbital::P);

                    coeff_j_i = coefficients(N,mo);

                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        atomicorbital=exponential;
                        coeff_j_i*=A;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+1,mo);

                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af;                       
                        Af = CalculateAngularNormalizationConstant(RenderOrbitals::PX);
                        CalculateAtomicOrbital(RenderOrbitals::PX,c,exponentialP,atomicorbital);
                        coeff_j_i=coeff_j_i*Ap*Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+2,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af;
                        Af = CalculateAngularNormalizationConstant(RenderOrbitals::PY);
                        CalculateAtomicOrbital(RenderOrbitals::PY,c,exponentialP,atomicorbital);
                        coeff_j_i=coeff_j_i*Ap*Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+3,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af;
                        Af = CalculateAngularNormalizationConstant(RenderOrbitals::PZ);
                        CalculateAtomicOrbital(RenderOrbitals::PZ,c,exponentialP,atomicorbital);
                        coeff_j_i=coeff_j_i*Ap*Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                }
                N+=4;
            }

            if (m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Type() == Orbital::P)
            {
                if ( (fabs(coefficients(N,mo)) > m_thresholdAO) || (fabs(coefficients(N+1,mo)) > m_thresholdAO) || (fabs(coefficients(N+2,mo)) > m_thresholdAO) )
                {
                    exponential.SetToZero();
                    for (size_t k=0; k<m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha().size(); ++k)
                    {
                        exponential.CalculateExponential(m_grid.L(),m_grid.L(),m_grid.L(),m_grid.Step(),c.x(),c.y(),c.z(),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha().at(k),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Xs().at(k));
                    }
                    float A = CalculateContractionConstant(m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Xs(),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha(),Orbital::P);
                    coeff_j_i = coefficients(N,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af;                        
                        Af = CalculateAngularNormalizationConstant(kryomol::RenderOrbitals::PX);
                        CalculateAtomicOrbital(RenderOrbitals::PX,c,exponential,atomicorbital);
                        coeff_j_i=coeff_j_i*A*Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+1,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af;                        
                        Af = CalculateAngularNormalizationConstant(RenderOrbitals::PY);
                        CalculateAtomicOrbital(RenderOrbitals::PY,c,exponential,atomicorbital);
                        coeff_j_i=coeff_j_i*A*Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+2,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af;                        
                        Af = CalculateAngularNormalizationConstant(RenderOrbitals::PZ);
                        CalculateAtomicOrbital(RenderOrbitals::PZ,c,exponential,atomicorbital);
                        coeff_j_i=coeff_j_i*A*Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                }
                N+=3;
            }


            if ( (m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Type() == Orbital::D) && (m_orbitaldata.TypeD()==6) )
            {
                if ( (fabs(coefficients(N,mo)) > m_thresholdAO) || (fabs(coefficients(N+1,mo)) > m_thresholdAO) || (fabs(coefficients(N+2,mo)) > m_thresholdAO) || (fabs(coefficients(N+3,mo)) > m_thresholdAO) || (fabs(coefficients(N+4,mo)) > m_thresholdAO) || (fabs(coefficients(N+5,mo)) > m_thresholdAO) )
                {
                    exponential.SetToZero();
                    for (size_t k=0; k<m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha().size(); ++k)
                    {
                        exponential.CalculateExponential(m_grid.L(),m_grid.L(),m_grid.L(),m_grid.Step(),c.x(),c.y(),c.z(),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha().at(k),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Xs().at(k));
                    }
                    float A = CalculateContractionConstant(m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Xs(),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha(),Orbital::D);
                    coeff_j_i = coefficients(N,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::DXX);
                        CalculateAtomicOrbital(RenderOrbitals::DXX,c,exponential,atomicorbital);
                        coeff_j_i=coeff_j_i*A*Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+1,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::DYY);
                        CalculateAtomicOrbital(RenderOrbitals::DYY,c,exponential,atomicorbital);
                        coeff_j_i=coeff_j_i*A*Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+2,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::DZZ);
                        CalculateAtomicOrbital(RenderOrbitals::DZZ,c,exponential,atomicorbital);
                        coeff_j_i=coeff_j_i*A*Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+3,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::DXY);
                        CalculateAtomicOrbital(RenderOrbitals::DXY,c,exponential,atomicorbital);
                        coeff_j_i=coeff_j_i*A*Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+4,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::DXZ);
                        CalculateAtomicOrbital(RenderOrbitals::DXZ,c,exponential,atomicorbital);
                        coeff_j_i=coeff_j_i*A*Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+5,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::DYZ);
                        CalculateAtomicOrbital(RenderOrbitals::DYZ,c,exponential,atomicorbital);
                        coeff_j_i=coeff_j_i*A*Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                }
                N+=6;
            }

            if ( (m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Type() == Orbital::D) && (m_orbitaldata.TypeD()==5) )
            {
                if ( (fabs(coefficients(N,mo)) > m_thresholdAO) || (fabs(coefficients(N+1,mo)) > m_thresholdAO) || (fabs(coefficients(N+2,mo)) > m_thresholdAO) || (fabs(coefficients(N+3,mo)) > m_thresholdAO) || (fabs(coefficients(N+4,mo)) > m_thresholdAO) )
                {
                    exponential.SetToZero();
                    for (size_t k=0; k<m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha().size(); ++k)
                    {
                        exponential.CalculateExponential(m_grid.L(),m_grid.L(),m_grid.L(),m_grid.Step(),c.x(),c.y(),c.z(),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha().at(k),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Xs().at(k));
                    }                    
                    float A = CalculateContractionConstant(m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Xs(),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha(),Orbital::D);
                    coeff_j_i = coefficients(N,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::DY0);
                        CalculateAtomicOrbital(RenderOrbitals::DY0,c,exponential,atomicorbital);
                        coeff_j_i*=A;
                        coeff_j_i*=Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+1,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::DY1);
                        CalculateAtomicOrbital(RenderOrbitals::DY1,c,exponential,atomicorbital);
                        coeff_j_i*=A;
                        coeff_j_i*=Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+2,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::DY2);
                        CalculateAtomicOrbital(RenderOrbitals::DY2,c,exponential,atomicorbital);
                        coeff_j_i*=A;
                        coeff_j_i*=Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+3,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::DY3);
                        CalculateAtomicOrbital(RenderOrbitals::DY3,c,exponential,atomicorbital);
                        coeff_j_i*=A;
                        coeff_j_i*=Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+4,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::DY4);
                        CalculateAtomicOrbital(RenderOrbitals::DY4,c,exponential,atomicorbital);
                        coeff_j_i*=A;
                        coeff_j_i*=Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                }
                N+=5;
            }

            if ( (m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Type() == Orbital::F) && (m_orbitaldata.TypeF()==10) )
            {
                if ( (fabs(coefficients(N,mo)) > m_thresholdAO) || (fabs(coefficients(N+1,mo)) > m_thresholdAO) || (fabs(coefficients(N+2,mo)) > m_thresholdAO) || (fabs(coefficients(N+3,mo)) > m_thresholdAO) || (fabs(coefficients(N+4,mo)) > m_thresholdAO) || (fabs(coefficients(N+5,mo)) > m_thresholdAO)  || (fabs(coefficients(N+6,mo)) > m_thresholdAO) || (fabs(coefficients(N+7,mo)) > m_thresholdAO) || (fabs(coefficients(N+8,mo)) > m_thresholdAO) || (fabs(coefficients(N+9,mo)) > m_thresholdAO))
                {
                    exponential.SetToZero();
                    for (size_t k=0; k<m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha().size(); ++k)
                    {
                        exponential.CalculateExponential(m_grid.L(),m_grid.L(),m_grid.L(),m_grid.Step(),c.x(),c.y(),c.z(),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha().at(k),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Xs().at(k));
                    }
                    float A = CalculateContractionConstant(m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Xs(),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha(),Orbital::F);
                    coeff_j_i = coefficients(N,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::FXXX);
                        CalculateAtomicOrbital(RenderOrbitals::FXXX,c,exponential,atomicorbital);
                        coeff_j_i=coeff_j_i*A*Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+1,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::FYYY);
                        CalculateAtomicOrbital(RenderOrbitals::FYYY,c,exponential,atomicorbital);
                        coeff_j_i=coeff_j_i*A*Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+2,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::FZZZ);
                        CalculateAtomicOrbital(RenderOrbitals::FZZZ,c,exponential,atomicorbital);
                        coeff_j_i=coeff_j_i*A*Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+3,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::FXYY);
                        CalculateAtomicOrbital(RenderOrbitals::FXYY,c,exponential,atomicorbital);
                        coeff_j_i=coeff_j_i*A*Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+4,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::FXXY);
                        CalculateAtomicOrbital(RenderOrbitals::FXXY,c,exponential,atomicorbital);
                        coeff_j_i=coeff_j_i*A*Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+5,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::FXXZ);
                        CalculateAtomicOrbital(RenderOrbitals::FXXZ,c,exponential,atomicorbital);
                        coeff_j_i=coeff_j_i*A*Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+6,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::FXZZ);
                        CalculateAtomicOrbital(RenderOrbitals::FXZZ,c,exponential,atomicorbital);
                        coeff_j_i=coeff_j_i*A*Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+7,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::FYZZ);
                        CalculateAtomicOrbital(RenderOrbitals::FYZZ,c,exponential,atomicorbital);
                        coeff_j_i=coeff_j_i*A*Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+8,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::FYYZ);
                        CalculateAtomicOrbital(RenderOrbitals::FYYZ,c,exponential,atomicorbital);
                        coeff_j_i=coeff_j_i*A*Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+9,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::FXYZ);
                        CalculateAtomicOrbital(RenderOrbitals::FXYZ,c,exponential,atomicorbital);
                        coeff_j_i=coeff_j_i*A*Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }

                }
                N+=10;
            }

            if ( (m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Type() == Orbital::F) && (m_orbitaldata.TypeF()==7) )
            {
                if ( (fabs(coefficients(N,mo)) > m_thresholdAO) || (fabs(coefficients(N+1,mo)) > m_thresholdAO) || (fabs(coefficients(N+2,mo)) > m_thresholdAO) || (fabs(coefficients(N+3,mo)) > m_thresholdAO) || (fabs(coefficients(N+4,mo)) > m_thresholdAO) || (fabs(coefficients(N+5,mo)) > m_thresholdAO)  || (fabs(coefficients(N+6,mo)) > m_thresholdAO) )
                {
                    exponential.SetToZero();
                    for (size_t k=0; k<m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha().size(); ++k)
                    {
                        exponential.CalculateExponential(m_grid.L(),m_grid.L(),m_grid.L(),m_grid.Step(),c.x(),c.y(),c.z(),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha().at(k),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Xs().at(k));
                    }
                    float A = CalculateContractionConstant(m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Xs(),m_orbitaldata.BasisCenters().at(j).Orbitals().at(i).Alpha(),Orbital::F);
                    coeff_j_i = coefficients(N,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::FY0);
                        CalculateAtomicOrbital(RenderOrbitals::FY0,c,exponential,atomicorbital);
                        coeff_j_i*=A;
                        coeff_j_i*=Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+1,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::FY1);
                        CalculateAtomicOrbital(RenderOrbitals::FY1,c,exponential,atomicorbital);
                        coeff_j_i*=A;
                        coeff_j_i*=Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+2,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::FY2);
                        CalculateAtomicOrbital(RenderOrbitals::FY2,c,exponential,atomicorbital);
                        coeff_j_i*=A;
                        coeff_j_i*=Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+3,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::FY3);
                        CalculateAtomicOrbital(RenderOrbitals::FY3,c,exponential,atomicorbital);
                        coeff_j_i*=A;
                        coeff_j_i*=Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+4,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::FY4);
                        CalculateAtomicOrbital(RenderOrbitals::FY4,c,exponential,atomicorbital);
                        coeff_j_i*=A;
                        coeff_j_i*=Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+5,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::FY5);
                        CalculateAtomicOrbital(RenderOrbitals::FY5,c,exponential,atomicorbital);
                        coeff_j_i*=A;
                        coeff_j_i*=Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                    coeff_j_i = coefficients(N+6,mo);
                    if (fabs(coeff_j_i) > m_thresholdAO)
                    {
                        float Af = CalculateAngularNormalizationConstant(RenderOrbitals::FY6);
                        CalculateAtomicOrbital(RenderOrbitals::FY6,c,exponential,atomicorbital);
                        coeff_j_i*=A;
                        coeff_j_i*=Af;
                        atomicorbital*=coeff_j_i;
                        AdditionAtomicOrbital(atom,atomicorbital,molecularorbital);
                    }
                }
                N+=7;
            }
        }
    }

    QApplication::restoreOverrideCursor();
#ifdef WITH_TIMERS
    //qDebug() << "Time elapsed: " << timer.elapsed() << endl;
    //std::cout << "Time elapsed: " << timer.elapsed() << std::endl;
#endif
}

void RenderOrbitals::AdditionAtomicOrbital(const Coordinate &distance, const OrbitalArray &atomicorbital, OrbitalArray &molecularorbital)
{
    //Update the total molecular orbital with the information of the molecular orbital around the atom
    int ix = floor((distance.x()-m_density.Origin().x())/m_grid.Step())-floor(m_grid.Nl()/2);
    int fx = floor((distance.x()-m_density.Origin().x())/m_grid.Step())+floor(m_grid.Nl()/2);
    int iy = floor((distance.y()-m_density.Origin().y())/m_grid.Step())-floor(m_grid.Nl()/2);
    int fy = floor((distance.y()-m_density.Origin().y())/m_grid.Step())+floor(m_grid.Nl()/2);
    int iz = floor((distance.z()-m_density.Origin().z())/m_grid.Step())-floor(m_grid.Nl()/2);
    int fz = floor((distance.z()-m_density.Origin().z())/m_grid.Step())+floor(m_grid.Nl()/2);
    //qDebug() << "CENTER SUBGRID GRID: " << ((floor((distance.x()-m_density.Origin().x())/m_grid.Step()))*m_grid.Step()+m_density.Origin().x()) << ((floor((distance.y()-m_density.Origin().y())/m_grid.Step()))*m_grid.Step()+m_density.Origin().y()) << ((floor((distance.z()-m_density.Origin().z())/m_grid.Step()))*m_grid.Step()+m_density.Origin().z()) << endl;

    int ax=0;
    for (int x=ix; x<fx;++x)
    {
        int ay=0;
        for (int y=iy; y<fy;++y)
        {
          int az=0;
          for (int z=iz; z<fz;++z)
          {
              if ( (x<(int)molecularorbital.NX()) && (y<(int)molecularorbital.NY()) && (z<(int)molecularorbital.NZ()) && (x>=0) && (y>=0) && (z>=0) )
              {
                  molecularorbital(x,y,z)=molecularorbital(x,y,z)+atomicorbital(ax,ay,az);
                  //qDebug() << "(" << x << "," << y << "," << z << ") -> " << "(" << ax << "," << ay << "," << az << ") -> " <<  atomicorbital(ax,ay,az) << molecularorbital(x,y,z) << endl;
              }
              ++az;
          }
          ++ay;
        }
        ++ax;
    }
}

float RenderOrbitals::CalculateContractionConstant(std::vector<float> &Xs, std::vector<float> &Alpha, Orbital::OrbitalType type)
{
    float N = 0.0;
    switch (type)
    {
    case Orbital::S:
        for (size_t i=0; i<Xs.size();++i)
            for (size_t j=0; j<Xs.size();++j)
                N += (Xs.at(i)*Xs.at(j))/(float)(sqrt((Alpha.at(i)+Alpha.at(j))*(Alpha.at(i)+Alpha.at(j))*(Alpha.at(i)+Alpha.at(j))));
        break;
    case Orbital::P:
        for (size_t i=0; i<Xs.size();++i)
            for (size_t j=0; j<Xs.size();++j)
                N += (Xs.at(i)*Xs.at(j))/(float)(sqrt((Alpha.at(i)+Alpha.at(j))*(Alpha.at(i)+Alpha.at(j))*(Alpha.at(i)+Alpha.at(j)))*(Alpha.at(i)+Alpha.at(j)));
        break;
    case Orbital::D:
        for (size_t i=0; i<Xs.size();++i)
            for (size_t j=0; j<Xs.size();++j)
                N += (Xs.at(i)*Xs.at(j))/(float)(sqrt((Alpha.at(i)+Alpha.at(j))*(Alpha.at(i)+Alpha.at(j))*(Alpha.at(i)+Alpha.at(j)))*(Alpha.at(i)+Alpha.at(j))*(Alpha.at(i)+Alpha.at(j)));
        break;
    case Orbital::F:
        for (size_t i=0; i<Xs.size();++i)
            for (size_t j=0; j<Xs.size();++j)
                N += (Xs.at(i)*Xs.at(j))/(float)(sqrt((Alpha.at(i)+Alpha.at(j))*(Alpha.at(i)+Alpha.at(j))*(Alpha.at(i)+Alpha.at(j)))*(Alpha.at(i)+Alpha.at(j))*(Alpha.at(i)+Alpha.at(j))*(Alpha.at(i)+Alpha.at(j)));
        break;
    default:
        N = 1.0;
        break;

    }

    N = (float) (sqrt(1/(sqrt((PC::Pi*PC::Pi*PC::Pi))*N)));

    return N;
}

float RenderOrbitals::CalculateAngularNormalizationConstant(Basis basis)
{
    float N;

    switch (basis)
    {
        case RenderOrbitals::PX:
            N = 0.5;
            break;
        case RenderOrbitals::PY:
            N = 0.5;
            break;
        case RenderOrbitals::PZ:
            N = 0.5;
            break;
        case RenderOrbitals::DXX:
            N = 0.75;
            break;
        case RenderOrbitals::DYY:
            N = 0.75;
            break;
        case RenderOrbitals::DZZ:
            N = 0.75;
            break;
        case RenderOrbitals::DXY:
            N = 0.25;
            break;
        case RenderOrbitals::DXZ:
            N = 0.25;
            break;
        case RenderOrbitals::DYZ:
            N = 0.25;
            break;
        case RenderOrbitals::FXXX:
            N = 1.875;
            break;
        case RenderOrbitals::FYYY:
            N = 1.875;
            break;
        case RenderOrbitals::FZZZ:
            N = 1.875;
            break;
        case RenderOrbitals::FXXY:
            N = 0.375;
            break;
        case RenderOrbitals::FXYY:
             N = 0.375;
            break;
        case RenderOrbitals::FXXZ:
            N = 0.375;
            break;
        case RenderOrbitals::FXZZ:
            N = 0.375;
            break;
        case RenderOrbitals::FYYZ:
            N = 0.375;
            break;
        case RenderOrbitals::FYZZ:
            N = 0.375;
            break;
        case RenderOrbitals::FXYZ:
            N = 0.125;
            break;
        case RenderOrbitals::DY0:
            N = 6.75;
            break;
        case RenderOrbitals::DY1:
            N = 0.25;
            break;
        case RenderOrbitals::DY2:
            N = 0.25;
            break;
        case RenderOrbitals::DY3:
            N = 2.25;
            break;
        case RenderOrbitals::DY4:
            N = 0.25;
            break;
        case RenderOrbitals::FY0:
            N = 16.875;
            break;
        case RenderOrbitals::FY1:
            N = 16.875;
            break;
        case RenderOrbitals::FY2:
            N = 16.875;
            break;
        case RenderOrbitals::FY3:
            N = 1.125;
            break;
        case RenderOrbitals::FY4:
            N = 0.125;
            break;
        case RenderOrbitals::FY5:
            N = 5.625;
            break;
        case RenderOrbitals::FY6:
            N = 5.625;
            break;
         default:
            N = 1.0;
            break;
    }

    N = (float) (sqrt(1/N));

    return N;
}


void RenderOrbitals::CalculateTotalDensity()
{
    #ifdef WITH_TIMERS
        QTime timerDensity;
        timerDensity.start();
    #endif

    QApplication::setOverrideCursor(QCursor(Qt::BusyCursor));

    m_totaldensity = OrbitalArray(m_density.Nx(),m_density.Ny(),m_density.Nz(),0);
    OrbitalArray molecularorbital(m_density.Nx(),m_density.Ny(),m_density.Nz());
    for (size_t i=0; i<m_orbitaldata.Coefficients().NColumns(); i++)
    {
        CalculateMolecularOrbital(i,molecularorbital,false);
        m_totaldensity.Hamard(molecularorbital,molecularorbital);
    }

    QApplication::restoreOverrideCursor();

#ifdef WITH_TIMERS
    //qDebug() << "Time elapsed for Total Density: " << timerDensity.elapsed() << endl;
    //std::cout << "Time elapsed for Total Density: " << timerDensity.elapsed() << std::endl;
#endif

}

void RenderOrbitals::CalculateSpinDensity()
{
    QApplication::setOverrideCursor(QCursor(Qt::BusyCursor));

    m_spindensity = OrbitalArray(m_density.Nx(),m_density.Ny(),m_density.Nz(),0);
    OrbitalArray betadensity(m_density.Nx(),m_density.Ny(),m_density.Nz(),0);
    OrbitalArray molecularorbital(m_density.Nx(),m_density.Ny(),m_density.Nz());
    OrbitalArray betaorbital(m_density.Nx(),m_density.Ny(),m_density.Nz());

    for (size_t i=0; i<m_orbitaldata.Coefficients().NColumns(); i++)
    {
        std::cout << "Molecular Orbital: " << i << " of " << m_orbitaldata.Coefficients().NColumns()-1 << std::endl;

        CalculateMolecularOrbital(i,molecularorbital,false);
        CalculateMolecularOrbital(i,betaorbital,true);

        m_spindensity.Hamard(molecularorbital,molecularorbital);
        betadensity.Hamard(betaorbital,betaorbital);
        qDebug() << molecularorbital(m_density.Nx()/2,m_density.Nx()/2,m_density.Nx()/2) << m_spindensity(m_density.Nx()/2,m_density.Nx()/2,m_density.Nx()/2) << endl;


    }

    m_spindensity-=betadensity;

    QApplication::restoreOverrideCursor();
}


void RenderOrbitals::CalculateHomo()
{
    QApplication::setOverrideCursor(QCursor(Qt::BusyCursor));

    m_homo = OrbitalArray(m_density.Nx(),m_density.Ny(),m_density.Nz());
    CalculateMolecularOrbital(m_orbitaldata.Homo()-1,m_homo,false);

    QApplication::restoreOverrideCursor();
}

void RenderOrbitals::CalculateLumo()
{
    QApplication::setOverrideCursor(QCursor(Qt::BusyCursor));

    m_lumo = OrbitalArray(m_density.Nx(),m_density.Ny(),m_density.Nz());
    if (m_beta)
        CalculateMolecularOrbital(m_orbitaldata.Lumo()-1,m_lumo,true);
    else
        CalculateMolecularOrbital(m_orbitaldata.Lumo()-1,m_lumo,false);

    QApplication::restoreOverrideCursor();
}


void RenderOrbitals::CalculateTransitionChange(size_t tc, OrbitalArray& transition)
{
    QApplication::setOverrideCursor(QCursor(Qt::BusyCursor));

    transition.SetToZero();

    OrbitalArray ground(m_density.Nx(),m_density.Ny(),m_density.Nz(),0);
    OrbitalArray excited(m_density.Nx(),m_density.Ny(),m_density.Nz(),0);

    std::vector<TransitionChange> transitiondata = m_transitiondata.at(tc);

    OrbitalArray orbital_i(m_density.Nx(),m_density.Ny(),m_density.Nz());;
    OrbitalArray orbital_j(m_density.Nx(),m_density.Ny(),m_density.Nz());;

    if (m_beta)
    {               
        for (size_t it=0; it<transitiondata.size(); ++it)
        {
            TransitionChange t = transitiondata.at(it);

            std::string& si = t.OrbitalSI();
            if (si.find("A")!=std::string::npos)
            {
                StringTokenizer token(si,"A");
                CalculateMolecularOrbital(atof(token.at(0))-1,orbital_i,false);
            }
            else
            {
                StringTokenizer token(si,"B");
                CalculateMolecularOrbital(atof(token.at(0))-1,orbital_i,true);
            }

            std::string& sj = t.OrbitalSJ();
            if (sj.find("A")!=std::string::npos)
            {
                StringTokenizer token(sj,"A");
                CalculateMolecularOrbital(atof(token.at(0))-1,orbital_j,false);
            }
            else
            {
                StringTokenizer token(sj,"B");
                CalculateMolecularOrbital(atof(token.at(0))-1,orbital_j,true);
            }

            orbital_i *= t.Coefficient();
            orbital_j *= t.Coefficient();

            ground += orbital_i;
            excited += orbital_j;
        }
    }
    else
    {
        for (size_t it=0; it<transitiondata.size(); ++it)
        {
            TransitionChange t = transitiondata.at(it);
            CalculateMolecularOrbital(t.OrbitalI()-1,orbital_i,false);
            CalculateMolecularOrbital(t.OrbitalJ()-1,orbital_j,false);

            orbital_i *= t.Coefficient();
            orbital_j *= t.Coefficient();

            ground += orbital_i;
            excited += orbital_j;
        }
    }

    transition.Hamard(ground,excited);

    QApplication::restoreOverrideCursor();
}

void RenderOrbitals::CalculateDensityChange(size_t tc, OrbitalArray& transition)
{
    QApplication::setOverrideCursor(QCursor(Qt::BusyCursor));

    transition.SetToZero();

    std::vector<TransitionChange> transitiondata = m_transitiondata.at(tc);

    OrbitalArray orbital_i(m_density.Nx(),m_density.Ny(),m_density.Nz());
    OrbitalArray orbital_j(m_density.Nx(),m_density.Ny(),m_density.Nz());

    if (m_beta)
    {
        for (size_t it=0; it<transitiondata.size(); ++it)
        {
            TransitionChange t = transitiondata.at(it);

            std::string& si = t.OrbitalSI();
            if (si.find("A")!=std::string::npos)
            {
                StringTokenizer token(si,"A");
                CalculateMolecularOrbital(atof(token.at(0))-1,orbital_i,false);
            }
            else
            {
                StringTokenizer token(si,"B");
                CalculateMolecularOrbital(atof(token.at(0))-1,orbital_i,true);
            }

            std::string& sj = t.OrbitalSJ();
            if (sj.find("A")!=std::string::npos)
            {
                StringTokenizer token(sj,"A");
                CalculateMolecularOrbital(atof(token.at(0))-1,orbital_j,false);
            }
            else
            {
                StringTokenizer token(sj,"B");
                CalculateMolecularOrbital(atof(token.at(0))-1,orbital_j,true);
            }
            transition.SquareDifference(t.Coefficient()*0.5,orbital_j,orbital_i);
        }
    }
    else
    {
        for (size_t it=0; it<transitiondata.size(); ++it)
        {
            TransitionChange t = transitiondata.at(it);
            CalculateMolecularOrbital(t.OrbitalI()-1,orbital_i,false);
            CalculateMolecularOrbital(t.OrbitalJ()-1,orbital_j,false);

            transition.SquareDifference(t.Coefficient()*0.5,orbital_j,orbital_i);
        }
    }

    QApplication::restoreOverrideCursor();
}


void RenderOrbitals::ShowHomo()
{
    if (m_homo.Empty())
        CalculateHomo();

    m_density.SetDensityMatrix(m_homo);
    m_density.RenderDensityData();
}

void RenderOrbitals::ShowLumo()
{
    if (m_lumo.Empty())
        CalculateLumo();

    m_density.SetDensityMatrix(m_lumo);
    m_density.RenderDensityData();
}

void RenderOrbitals::ShowTotalDensity()
{
    if (m_totaldensity.Empty())
        CalculateTotalDensity();

    m_density.SetDensityMatrix(m_totaldensity);
    m_density.RenderDensityData();
}

void RenderOrbitals::ShowSpinDensity()
{
    if (m_spindensity.Empty())
        CalculateSpinDensity();

    m_density.SetDensityMatrix(m_spindensity);
    m_density.RenderDensityData();
}

void RenderOrbitals::ShowMolecularOrbital(size_t mo)
{
    OrbitalArray array(m_density.Nx(),m_density.Ny(),m_density.Nz());

    CalculateMolecularOrbital(mo,array,false);

    m_density.SetDensityMatrix(array);
    m_density.RenderDensityData();
}

void RenderOrbitals::ShowBetaMolecularOrbital(size_t mo)
{
    OrbitalArray array(m_density.Nx(),m_density.Ny(),m_density.Nz());

    CalculateMolecularOrbital(mo,array,true);

    m_density.SetDensityMatrix(array);
    m_density.RenderDensityData();
}

void RenderOrbitals::ShowTransitionChange(size_t tc)
{
    OrbitalArray array(m_density.Nx(),m_density.Ny(),m_density.Nz());

    CalculateTransitionChange(tc,array);

    m_density.SetDensityMatrix(array);
    m_density.RenderDensityData();
}

void RenderOrbitals::ShowDensityChange(size_t tc)
{
    OrbitalArray array(m_density.Nx(),m_density.Ny(),m_density.Nz());

    CalculateDensityChange(tc,array);

    m_density.SetDensityMatrix(array);
    m_density.RenderDensityData();
}

void RenderOrbitals::CleanMatrices(bool b_homo, bool b_lumo, bool b_total, bool b_spin, int b_orbital, int b_betaorbital)
{
    m_homo.Clear();
    m_lumo.Clear();
    m_totaldensity.Clear();
    m_spindensity.Clear();

    if (b_homo)
    {
        CalculateHomo();

        m_density.SetDensityMatrix(m_homo);
        m_density.RenderDensityData();
    }
    if (b_lumo)
    {
        CalculateLumo();

        m_density.SetDensityMatrix(m_lumo);
        m_density.RenderDensityData();
    }
    if (b_total)
    {
        CalculateTotalDensity();

        m_density.SetDensityMatrix(m_totaldensity);
        m_density.RenderDensityData();
    }
    if ((b_spin)&&(m_beta))
    {
        CalculateSpinDensity();

        m_density.SetDensityMatrix(m_spindensity);
        m_density.RenderDensityData();
    }
    if (b_orbital)
    {
        OrbitalArray array(m_density.Nx(),m_density.Ny(),m_density.Nz());

        CalculateMolecularOrbital(b_orbital,array,false);

        m_density.SetDensityMatrix(array);
        m_density.RenderDensityData();
    }
    if ((b_betaorbital)&&(m_beta))
    {
        OrbitalArray array(m_density.Nx(),m_density.Ny(),m_density.Nz());

        CalculateMolecularOrbital(b_betaorbital,array,true);

        m_density.SetDensityMatrix(array);
        m_density.RenderDensityData();
    }
}
