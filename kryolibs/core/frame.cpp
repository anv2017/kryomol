/*****************************************************************************************
                            frame.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <limits>


#include "frame.h"
#include "molecule.h"
#include "ringperceptor.h"
#include "grid.h"


struct fcolor {
    float h;
    float s;
    float l;
};

class kryomol::FramePrivate
{
public:
    FramePrivate() {}
    ~FramePrivate() {}
    std::vector<Bond> m_bonds;
    std::vector<Coordinate> m_xyz;
    std::vector<Coordinate> m_gradient;
    Energy m_kenergy;
    Energy m_venergy;
    Energy m_rmsforce;
    double m_maximumforce;
    double m_maximumdisplacement;
    double m_rmsdisplacement;
    double m_S2;

    Threshold m_threshold;
    Coordinate m_dipole;
    OrbitalData m_orbitaldata;
    ElectronicDensity m_electronicdensity;

    Grid m_grid;

    D1Array<double> m_forces;
    D2Array<double> m_hessian;
    D2Array<double> m_dipolederivatives;
    std::vector< D2Array<double> > m_cshifttensors;
    std::vector<double> m_espcharges;
    std::vector<Frequency>  m_modes;
    std::vector<Coordinate> m_inputorientation;
    std::vector<Spectralline> m_spectrallines;
    std::vector<RenderDensity> m_positivedensity;
    std::vector<RenderDensity> m_negativedensity;
    std::vector< std::vector<TransitionChange> > m_transitionchanges;
    fcolor m_color;
};

using namespace kryomol;

Frame::Frame ( Molecule* molecule ) : m_molecule ( molecule )
{
    m_private = new FramePrivate();
}

Frame::~Frame()
{
    delete m_private;
}

Frame::Frame ( const Frame& frame )
{
    if ( this != &frame )
    {
        m_private = new FramePrivate ( * ( frame.m_private ) );
    }
    m_molecule=frame.m_molecule;
}

Frame& Frame::operator= ( const Frame& frame )
{
    delete m_private;
    m_private = new FramePrivate ( * ( frame.m_private ) );
    m_molecule=frame.m_molecule;
    return  *this;
}

Coordinate Frame::MassCenter() const
{
    std::vector<Atom>::const_iterator it;
    std::vector<Coordinate>::const_iterator ct;
    Coordinate cm;
    double total=0;
    for ( it=m_molecule->Atoms().begin(),ct=m_private->m_xyz.begin();it!=m_molecule->Atoms().end();++it,++ct )
    {
        cm.x() +=ct->x() *it->AtomicMass();
        cm.y() +=ct->y() *it->AtomicMass();
        cm.z() +=ct->z() *it->AtomicMass();
        total+=it->AtomicMass();
    }

    cm/=total;
    return cm;
}


D2Array<double> Frame::InertiaTensor() const
{

    D2Array<double> imatrix ( 3,3,0. );

    std::vector<Coordinate>::const_iterator at;
    std::vector<Coordinate> cmcoords;
    cmcoords.resize ( m_private->m_xyz.size() );
    std::vector<Coordinate>::iterator ct;

    Coordinate mc=MassCenter();
    for ( at=m_private->m_xyz.begin(),ct=cmcoords.begin();at!=m_private->m_xyz.end();++at,++ct )
    {
        ( *ct ) = ( *at )-mc;
    }
    std::vector<Atom>::const_iterator mat;
    for ( mat=m_molecule->Atoms().begin(),ct=cmcoords.begin();mat!=m_molecule->Atoms().end();++mat,++ct )
    {
        imatrix ( 0,0 ) +=mat->AtomicMass() * ( ct->y() *ct->y() +ct->z() *ct->z() );
        imatrix ( 1,1 ) +=mat->AtomicMass() * ( ct->x() *ct->x() +ct->z() *ct->z() );
        imatrix ( 2,2 ) +=mat->AtomicMass() * ( ct->x() *ct->x() +ct->y() *ct->y() );

        imatrix ( 0,1 ) =imatrix ( 1,0 )-=mat->AtomicMass() * ( ct->x() *ct->y() );
        imatrix ( 0,2 ) =imatrix ( 2,0 )-=mat->AtomicMass() * ( ct->x() *ct->z() );

        imatrix ( 1,2 ) =imatrix ( 2,1 )-=mat->AtomicMass() * ( ct->y() *ct->z() );

    }

    return imatrix;
}

D2Array<double> Frame::GyrationTensor() const
{
    D2Array<double> imatrix ( 3,3,0. );

    std::vector<Coordinate>::const_iterator at;
    std::vector<Coordinate> cmcoords;
    cmcoords.resize ( m_private->m_xyz.size() );
    std::vector<Coordinate>::iterator ct;

    Coordinate mc=Centroid();
    for ( at=m_private->m_xyz.begin(),ct=cmcoords.begin();at!=m_private->m_xyz.end();++at,++ct )
    {
        ( *ct ) = ( *at )-mc;
    }
    std::vector<Atom>::const_iterator mat;
    for ( mat=m_molecule->Atoms().begin(),ct=cmcoords.begin();mat!=m_molecule->Atoms().end();++mat,++ct )
    {

        imatrix ( 0,0 ) +=ct->x() *ct->x();
        imatrix ( 1,1 ) +=ct->y() *ct->y();
        imatrix ( 2,2 ) +=ct->z() *ct->z();

        imatrix ( 0,1 ) =imatrix ( 1,0 ) += ( ct->x() *ct->y() );
        imatrix ( 0,2 ) =imatrix ( 2,0 ) += ( ct->x() *ct->z() );

        imatrix ( 1,2 ) =imatrix ( 2,1 ) += ( ct->y() *ct->z() );


    }
    for ( int i=0;i<3;++i )
        for ( int j=0;j<3;++j )
            imatrix ( i,j ) /=m_molecule->Atoms().size();

    return imatrix;
}
float Frame::Angle ( int a1, int a2, int a3 ) const
{

    Coordinate v1=m_private->m_xyz[a1] - m_private->m_xyz[a2];
    Coordinate v2=m_private->m_xyz[a3] - m_private->m_xyz[a2];

    float angle=Coordinate::Angle ( v1,v2 );
    return angle;
}

std::vector < std::vector<size_t> > Frame::Rings ( size_t size/*=0*/ ) const
{
    RingPerceptor pr ( this );
    return pr.Rings ( size );
}

Coordinate Frame::Centroid() const
{
    std::vector<Coordinate>::const_iterator it;
    Coordinate c;
    for ( it=m_private->m_xyz.begin();it!=m_private->m_xyz.end();++it )
    {
        c+= ( *it );
    }
    c/=m_private->m_xyz.size();
    return c;
}

/** Get the a->b vector */
Coordinate Frame::Vector ( int i,int j ) const
{

    return m_private->m_xyz[i]-m_private->m_xyz[j];
}

float Frame::Dihedral ( int i, int j, int k, int l ) const
{
    return Coordinate::Dihedral ( m_private->m_xyz[i],m_private->m_xyz[j],m_private->m_xyz[k],m_private->m_xyz[l] );
}

const std::vector<Bond>& Frame::Bonds() const
{
    return m_private->m_bonds;
}
std::vector<Bond>& Frame::Bonds()
{
    return m_private->m_bonds;
}

const std::vector<Coordinate>& Frame::XYZ() const
{
    return m_private->m_xyz;
}

std::vector<Coordinate>& Frame::XYZ()
{
    return m_private->m_xyz;
}

const std::vector<Coordinate>& Frame::Gradient() const
{
    return m_private->m_gradient;
}

std::vector<Coordinate>& Frame::Gradient()
{
    return m_private->m_gradient;
}

Energy& Frame::KineticEnergy()
{
    return m_private->m_kenergy;

}

const Energy& Frame::KineticEnergy() const
{
    return m_private->m_kenergy;
}

const Energy& Frame::PotentialEnergy() const
{
    return m_private->m_venergy;
}

Energy& Frame::PotentialEnergy()
{
    return m_private->m_venergy;
}

Energy Frame::TotalEnergy() const
{
    return m_private->m_venergy + m_private->m_kenergy;
}

const Energy& Frame::RMSForce() const
{
    return m_private->m_rmsforce;
}

Energy& Frame::RMSForce()
{
    return m_private->m_rmsforce;
}

float Frame::Distance ( int i, int j ) const
{
    return Coordinate::Distance ( m_private->m_xyz[i],m_private->m_xyz[j] );
}

void Frame::SetDihedral ( size_t i,size_t j, size_t k, size_t l, float dihedral )
{
    float oldd= Dihedral ( i,j,k,l );
    float pass=dihedral-oldd;

    //m_rotaxis[0]=m_atom[j].c;
    //m_rotaxis[1]=m_atom[k].c;
    RotateBond ( j,k,pass,1 );
}

void Frame::RotateBond ( size_t i,size_t j,float angle, bool clockwise )
{
    std::vector<size_t> rotatedatoms;
    RotateNeighbours ( XYZ() [i],XYZ() [j],i,j,angle,clockwise,rotatedatoms );

}

void Frame::RotateNeighbours ( const Coordinate& axisorigin,const Coordinate& axisend,size_t i,size_t j,float pass,bool clockwise,std::vector<size_t>& rotatedatoms )
{
    Coordinate& atom=XYZ() [j];

    std::vector<size_t>::iterator pos=std::find ( rotatedatoms.begin(), rotatedatoms.end() , j );
    if ( pos == rotatedatoms.end() )
    {

        if (clockwise )
        {
            atom=Coordinate::RotAroundAxis ( atom,axisorigin,axisend,pass );
        }
        else
        {
            atom=Coordinate::RotAroundAxis ( atom,axisorigin,axisend,-pass );
        }
        rotatedatoms.push_back ( j );
    }
    else return;


    const std::vector<Bond>& bonds=ParentMolecule()->Bonds();
    std::vector<Bond>::const_iterator it;
    for ( it=bonds.begin();it!=bonds.end();++it )
    {

        if ( it->I() == j && it->J() != i )
            RotateNeighbours ( axisorigin,axisend,j, it->J(), pass,clockwise,rotatedatoms );
        if ( it->J() == j && it->I() != i )
            RotateNeighbours ( axisorigin,axisend,j, it->I(), pass,clockwise,rotatedatoms );
    }
}

std::vector< D2Array<double> >&  Frame::CShiftTensors()
{
    return m_private->m_cshifttensors;
}

const std::vector< D2Array<double> >&  Frame::CShiftTensors() const
{
    return m_private->m_cshifttensors;
}

OrbitalData& Frame::OrbitalsData()
{
    return m_private->m_orbitaldata;
}

const OrbitalData& Frame::OrbitalsData() const
{
    return m_private->m_orbitaldata;
}

std::vector< std::vector<TransitionChange> >& Frame::TransitionChanges()
{
    return m_private->m_transitionchanges;
}

const std::vector< std::vector<TransitionChange> >& Frame::TransitionChanges() const
{
    return m_private->m_transitionchanges;
}

ElectronicDensity& Frame::ElectronicDensityData()
{
    return m_private->m_electronicdensity;
}

const ElectronicDensity& Frame::ElectronicDensityData() const
{
    return m_private->m_electronicdensity;
}

std::vector<RenderDensity>& Frame::PositiveDensity()
{
    return m_private->m_positivedensity;
}

const std::vector<RenderDensity>& Frame::PositiveDensity() const
{
    return m_private->m_positivedensity;
}

std::vector<RenderDensity>& Frame::NegativeDensity()
{
    return m_private->m_negativedensity;
}

const std::vector<RenderDensity>& Frame::NegativeDensity() const
{
    return m_private->m_negativedensity;
}

Grid& Frame::GetGrid()
{
    return m_private->m_grid;
}

const Grid& Frame::GetGrid() const
{
    return m_private->m_grid;
}

std::ostream& kryomol::operator << ( std::ostream& s,const kryomol::Frame& m )
{
    std::vector<Atom>::const_iterator it;
    std::vector<Coordinate>::const_iterator ct;
    for ( it=m.ParentMolecule()->Atoms().begin(),ct=m.XYZ().begin();it!=m.ParentMolecule()->Atoms().end();++it,++ct )
    {

        s << std::setiosflags ( std::ios::fixed ) << std::setw ( 3 ) << std::setiosflags ( std::ios::left ) << it->Symbol()
          << std::setw ( 12 )  << std::resetiosflags ( std::ios::left )
          << std::setiosflags ( std::ios::right ) << std::setiosflags ( std::ios::fixed ) << std::setprecision ( 6 )
          << ct->x() <<
             std::setw ( 12 ) << ct->y()
          << std::setw ( 12 ) << ct->z()
          << std::resetiosflags ( std::ios::right ) << std::endl;
    }

    s << std::resetiosflags ( std::ios::fixed );
    return s;
}



std::vector<double>* Frame::Charges(Frame::Charge type)
{
    switch( type )
    {
    case ESP:
        return &m_private->m_espcharges;
    default:
        return NULL;
    }
    return NULL;
}

const std::vector<double>* Frame::Charges(Frame::Charge type) const
{
    switch( type )
    {
    case ESP:
        return &m_private->m_espcharges;
    default:
        return NULL;
    }
    return NULL;
}

void Frame::AllocateVectors(bool rottransmodes/*=false*/)
{
    size_t size;
    rottransmodes ? size =3*ParentMolecule()->Atoms().size() : size = 3*ParentMolecule()->Atoms().size()-6;

    m_private->m_modes.reserve(size);
    ParentMolecule()->GetMode().resize(size);

    for( std::vector< std::vector<Coordinate> >::iterator it=ParentMolecule()->GetMode().begin();it != ParentMolecule()->GetMode().end();++it)
    {
        it->resize(ParentMolecule()->Atoms().size());
    }

}


void Frame::AllocateForces()
{

    m_private->m_forces.Initialize(3*ParentMolecule()->Atoms().size());
}

void Frame::AllocateHessian()
{

    m_private->m_hessian.Initialize(3*ParentMolecule()->Atoms().size(),3*ParentMolecule()->Atoms().size());
    m_private->m_dipolederivatives.Initialize(3*ParentMolecule()->Atoms().size(),3);

}

void Frame::SetEnergy(double e, const std::string& level)
{
    m_private->m_venergy=e; ParentMolecule()->SetEnergyLevel(level);
}

double Frame::GetEnergy() const
{
    return m_private->m_venergy.Value();
}


std::pair<Coordinate,Coordinate> Frame::Box() const
{

    std::pair<Coordinate,Coordinate> box;
    /*box.first.x=std::numeric_limits<double>::max(); // lowest x
    box.second.x=std::numeric_limits<double>::min(); //highest x
    box.first.y=std::numeric_limits<double>::max(); //lowest y
    box.second.y=std::numeric_limits<double>::min();  //highest t
    box.first.z=std::numeric_limits<double>::max(); //lowest z
    box.second.z=std::numeric_limits<double>::min();  //highest z
    for (std::vector<Coordinate>::const_iterator it=m_private->m_xyz.begin();it!=m_private->m_xyz.end();++it)
    {
        if ( box.first.x > it->x ) box.first.x=it->x;
        if ( box.second.x < it->x ) box.second.x=it->x;
        if ( box.first.y > it->y ) box.first.y = it->y;
        if ( box.second.y < it->y ) box.second.y=it->y;
        if ( box.first.z > it->z ) box.first.z=it->z;
        if ( box.second.z < it->z ) box.second.z=it->z;

    }*/

    return box;

}

void Frame::CalculateGrid(float step)
{
    Coordinate centroid = this->Centroid();

    float xmax=0; float ymax=0; float zmax=0;
    for (size_t i=0; i<m_private->m_xyz.size(); i++)
    {
        if (xmax < abs(m_private->m_xyz.at(i).x()-centroid.x()))
            xmax = abs(m_private->m_xyz.at(i).x()-centroid.x());
        if (ymax < abs(m_private->m_xyz.at(i).y()-centroid.y()))
            ymax = abs(m_private->m_xyz.at(i).y()-centroid.y());
        if (zmax < abs(m_private->m_xyz.at(i).z()-centroid.z()))
            zmax = abs(m_private->m_xyz.at(i).z()-centroid.z());
    }

    // Minimum size of the grid: 2 Angstrom
    if (xmax<2)  xmax=2;
    if (ymax<2)  ymax=2;
    if (zmax<2)  zmax=2;

    //Both sides of the grid
    xmax = 2*xmax;
    ymax = 2*ymax;
    zmax = 2*zmax;

    //Increase the size of the grid for representing small values of isovalue
    xmax = 4*xmax;
    ymax = 4*ymax;
    zmax = 4*zmax;

    //The size of the default subgrid is 6.0 Angstrom
    float lmax = 6.0;

    qDebug() << "Grid dimension: x=" << xmax << " y=" << ymax << " z=" << zmax << " l=" << lmax <<endl;

    m_private->m_grid = Grid(xmax,ymax,zmax,step,lmax);
}

void Frame::SetColor(float h, float s, float l)
{
    m_private->m_color.h=h;
    m_private->m_color.s=s;
    m_private->m_color.l=l;
}

void Frame::GetColor(float& h, float& s, float& l) const
{
    h=m_private->m_color.h;
    s=m_private->m_color.s;
    l=m_private->m_color.l;
}


void Frame::SetS2(double e) { m_private->m_S2=e; }
void Frame::SetRMSForce(double f) { m_private->m_rmsforce=f; }
void Frame::SetMaximumForce(double f) { m_private->m_maximumforce=f; }
void Frame::SetRMSDisplacement(double d) { m_private->m_rmsdisplacement=d; }
void Frame::SetMaximumDisplacement(double d) { m_private->m_maximumdisplacement=d; }
void Frame::SetThreshold(const Threshold& thr) { m_private->m_threshold=thr; }
void Frame::SetDipole(const Coordinate& coor) { m_private->m_dipole=coor; }
void Frame::SetGrid(const Grid &grid) { m_private->m_grid=grid; }
void Frame::SetOrbitalData(const OrbitalData &orbitaldata) { m_private->m_orbitaldata=orbitaldata; }
void Frame::SetTransitionChanges(const std::vector< std::vector<TransitionChange> > &transitions) { m_private->m_transitionchanges=transitions; }
void Frame::SetElectronicDensityData(const ElectronicDensity &density) {m_private->m_electronicdensity=density; }
void Frame::SetPositiveDensity(const std::vector<RenderDensity> &positivedensity) {m_private->m_positivedensity=positivedensity;}
void Frame::SetNegativeDensity(const std::vector<RenderDensity> &negativedensity) {m_private->m_negativedensity=negativedensity;}

double Frame::GetS2() const { return m_private->m_S2; }
double Frame::GetRMSForce() const { return m_private->m_rmsforce.Value(); }
double Frame::GetMaximumForce() const { return m_private->m_maximumforce; }
double Frame::GetRMSDisplacement() const { return m_private->m_rmsdisplacement; }
double Frame::GetMaximumDisplacement() const { return m_private->m_maximumdisplacement; }
Threshold Frame::GetThreshold() const { return m_private->m_threshold; }
Coordinate Frame::GetDipole() const { return m_private->m_dipole; }

D1Array<double>& Frame::GetForces()  { return m_private->m_forces; }
D2Array<double>& Frame::GetHessian()  { return m_private->m_hessian; }
std::vector<Frequency> &Frame::GetFrequencies() { return m_private->m_modes; }
const std::vector<Frequency>& Frame::GetFrequencies() const { return m_private->m_modes; }
std::vector<Spectralline>& Frame::GetSpectralLines() { return m_private->m_spectrallines; }
const std::vector<Spectralline>& Frame::GetSpectralLines() const { return m_private->m_spectrallines; }
