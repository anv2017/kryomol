/*****************************************************************************************
                            molecule.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <math.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "molecule.h"
#include "ringperceptor.h"
#include "stringtools.h"
#include "exception.h"


#ifndef M_PI
#define M_PI       3.14159265358979323846  // matches value in gcc v2 math.h
#endif


using namespace kryomol;

class kryomol::MoleculePrivate
{
public:

    MoleculePrivate()  : m_currentframe ( 0 )
    {}
    MoleculePrivate ( const MoleculePrivate& mol )
    {
        m_atoms=mol.m_atoms;
        m_currentframe=mol.m_currentframe;
        m_bonds=mol.m_bonds;
        m_frames=mol.m_frames;
        m_populations=mol.m_populations;
        m_residues.reserve ( mol.m_residues.size() );
        for ( std::vector<PDBResidue*>::const_iterator it=mol.m_residues.begin();it!=mol.m_residues.end();++it )
        {
            PDBResidue* res=new PDBResidue ( * ( *it ) );
            m_residues.push_back ( res );
        }
        m_energylevel=mol.m_energylevel;
    }
    MoleculePrivate& operator= ( const MoleculePrivate& mol )
    {
        if ( this != &mol )
        {
            m_atoms=mol.m_atoms;
            m_currentframe=mol.m_currentframe;
            m_bonds=mol.m_bonds;
            m_frames=mol.m_frames;
            m_populations=mol.m_populations;
            for ( std::vector<PDBResidue*>::iterator it=m_residues.begin();it!=m_residues.end();++it )
            {
                delete ( *it );
            }
            m_residues.clear();
            for ( std::vector<PDBResidue*>::const_iterator it=mol.m_residues.begin();it!=mol.m_residues.end();++it )
            {
                PDBResidue* res=new PDBResidue ( * ( *it ) );
                m_residues.push_back ( res );
            }
            m_energylevel = mol.m_energylevel;
        }
        return *this;
    }
    ~MoleculePrivate()
    {
        for ( std::vector<PDBResidue*>::iterator it=m_residues.begin();it!=m_residues.end();++it )
            delete ( *it );
    }
    std::vector<Atom> m_atoms;
    size_t m_currentframe;
    std::vector<Bond> m_bonds;
    std::vector<Frame> m_frames;
    std::vector<double> m_populations;
    std::vector<PDBResidue*> m_residues;
    std::string m_energylevel;
};


Molecule::Molecule()
{
    m_private = new MoleculePrivate();
}

Molecule::Molecule ( const Molecule& mol )
{
    m_private = new MoleculePrivate ( * ( mol.m_private ) );

}

Molecule& Molecule::operator= ( const Molecule& mol )
{
    if ( this != &mol )
    {
        delete m_private;
        m_private = new MoleculePrivate ( * ( mol.m_private ) );
    }
    return  *this;

}

Molecule::~Molecule()
{
    delete m_private;
}

Molecule Molecule::At ( size_t frame ) const
{
    Molecule molecule;
    molecule.Atoms() =this->Atoms();
    molecule.Bonds() =this->Bonds();
    molecule.Frames().push_back ( this->Frames().at ( frame ) );
    return molecule;
}


/** No descriptions */
// bool CMolecule::AreSp3Connected(int i, int j)
// {
//
//
//   Atom& H1=m_atoms[i];
//   Atom& H2=m_atoms[j];
//   Atom& C1=this->m_atoms[H1.conexion[0].index];
//   Atom& C2=this->m_atoms[H2.conexion[0].index];
//
//   if (C1.conexion.size()!=4 || C2.conexion.size()!=4)
//   {
//     std::cout << "No sp3 atoms";
//     return false;
//   }
//
//   bool b=false;
//   for (int k=0;k<4;k++)
//   {
//
//     if(C1.conexion[k].index==H2.conexion[0].index)
//       b=true;
//   }
//
//   return b;
// }


/** Calculate Molecular Weight */
double Molecule::CalculateWeight() const
{

    std::vector<Atom>::const_iterator it;
    float sum=0.0f;
    for ( it=m_private->m_atoms.begin();it!=m_private->m_atoms.end();++it )
    {
        sum+=it->AtomicMass();
    }

    return sum;

}

/** Calculate an approximate center of masses for better rotation handling
 it is not really CM cause all atoms have the same weight */
void Molecule::MoveToCentroid()
{


    std::vector<Frame>::iterator ft;
    for ( ft=m_private->m_frames.begin();ft!=m_private->m_frames.end();++ft )
    {
        Coordinate cm=ft->Centroid();
        std::vector<Coordinate>::iterator ct;
        for ( ct=ft->XYZ().begin();ct!=ft->XYZ().end();++ct )
        {
            ct->x()-=cm.x();
            ct->y()-=cm.y();
            ct->z()-=cm.z();
        }
    }


}

std::vector < std::vector<size_t> > Molecule::Rings ( size_t size/*=0*/ ) const
{
    std::cout << "Rings" << std::endl;
    RingPerceptor pr ( this );
    return pr.Rings ( size );
}


std::vector<int> Molecule::FindAtomByName ( const std::string& name ) const
{

    std::vector<Atom>::const_iterator it;

    int counter=0;
    std::vector<int> list;
    list.reserve ( m_private->m_atoms.size() );
    for ( it=m_private->m_atoms.begin();it!=m_private->m_atoms.end();++it )
    {
        if ( it->Symbol() == name )
            list.push_back ( counter );
        counter++;
    }

    return list;
}

/** Get Hybridization type for atom i */
Atom::hybridization Molecule::Hybridization ( size_t atom ) const
{
    if ( m_private->m_atoms.at ( atom ).Z() == 6 )
    {
        std::vector<size_t> ng=Neighbours ( atom );
        switch ( ng.size() )
        {
        case 4:
            return Atom::SP3;
        case 3:
            return Atom::SP2;
        case 2:
        case 1:
            return Atom::SP;
        default:
            return Atom::UNKNOWN;
        }
    }
    else
    {
        std::cout << "Not implemented for atoms other than C";
        return Atom::UNKNOWN;
    }

}

void Molecule::SetBonds ( bool eachframe/*=false*/ )
{

    m_private->m_bonds.clear();
    std::vector<Frame>::iterator ft;
    for ( ft=m_private->m_frames.begin();ft!=m_private->m_frames.end();++ft )
    {
        ft->Bonds().clear();
    }

    std::vector<Atom>::const_iterator at,at1;
    std::vector<Coordinate>::iterator ct,ct1;

    int counter1,counter2;
    std::vector<Bond>* conexion;
    if ( eachframe )
        ft=Frames().begin();
    else
    {
        conexion=&Bonds();
        ft=Frames().end()-1;
    }
    conexion->reserve ( 2*Atoms().size() ); //a hydrocarbon have nearly 2*N bonds
    for ( ;ft!=Frames().end();++ft )
    {
        if ( eachframe ) conexion=&ft->Bonds();
        for ( at=Atoms().begin(),counter1=0,ct=ft->XYZ().begin();at!=Atoms().end();++at,++ct,++counter1 )
        {
            for ( at1=at+1,counter2=counter1+1,ct1=ct+1;at1!=Atoms().end();++at1,++ct1,++counter2 )
            {
                if ( at->Z() > 0 && at1->Z() > 0 )
                {
                    if ( ( ct->z() - ct1->z() )  < 3.0 && ( ct->y() - ct1->y() )  < 3.0 && ( ct->x() - ct1->x() )  < 3.0 )
                    {
                        float d=Coordinate::Distance ( *ct,*ct1 );
                        if ( d < 3.0f )
                        {
                            //   if ( at->CovalentRadius()+at1->CovalentRadius() < d*1.1 )
                            //    conexion->push_back(Bond(counter1,counter2,Bond::SINGLE));
                            if ( at->Z() == 1 || at1->Z() == 1 )
                            {
                                if ( d <1.5f )
                                {
                                    conexion->push_back ( Bond ( counter1,counter2,Bond::SINGLE ) );
                                    //    std::cout << "connected(" <<counter1+1 << "," << counter2+1 << ")" << std::endl;
                                }
                            }

                            else if ( at->Z() <11 && at1->Z() < 11 )
                            {
                                if ( d< 2.0f )
                                {
                                    conexion->push_back ( Bond ( counter1,counter2,Bond::SINGLE ) );
                                    //   std::cout << "connected(" <<counter1+1 << "," << counter2+1 << ")" << std::endl;
                                }
                            }

                            else if ( ( at->Z() <11 && at1->Z() >= 11 ) ||
                                      ( at1->Z() <11 && at->Z() >= 11 ) )
                            {
                                if ( d < 2.3f )
                                    conexion->push_back ( Bond ( counter1,counter2,Bond::SINGLE ) );
                                //std::cout << "connected(" <<counter1+1 << "," << counter2+1 << ")" << std::endl;
                            }

                            else
                            {
                                conexion->push_back ( Bond ( counter1,counter2,Bond::SINGLE ) );
                            }
                        }
                    }
                }
            }
        }

    }
}


std::ostream& operator << ( std::ostream& s,const Molecule& m )
{   
    std::vector<Frame>::const_iterator ft;
    size_t n=0;
    for ( ft=m.Frames().begin();ft!=m.Frames().end();++ft )
    {
        s << "Frame # " << ++n << std::endl;
        s << ( *ft );
        s << std::endl;
    }

    return s;
}


std::vector<size_t> Molecule::Neighbours ( size_t i ) const
{
    std::vector<size_t> n;
    std::vector<Bond>::const_iterator bt;
    const std::vector<Bond>* b;
    if ( Bonds().empty() ) b=&CurrentFrame().Bonds();
    else b=&Bonds();
    for ( bt=b->begin();bt!=b->end();++bt )
    {
        if ( bt->I() == i ) n.push_back ( bt->J() );
        if ( bt->J() == i ) n.push_back ( bt->I() );
    }
    return n;
}

const std::vector<Atom>& Molecule::Atoms() const
{
    return m_private->m_atoms;
}
std::vector<Atom>& Molecule::Atoms()
{
    return m_private->m_atoms;
}

const std::vector<Frame>& Molecule::Frames() const
{
    return m_private->m_frames;
}

std::vector<Frame>& Molecule::Frames()
{
    return m_private->m_frames;
}

const Frame& Molecule::CurrentFrame() const
{
    return m_private->m_frames.at ( m_private->m_currentframe );
}

Frame& Molecule::CurrentFrame()
{
    return m_private->m_frames.at ( m_private->m_currentframe );
}

size_t Molecule::CurrentFrameIndex() const
{
    return m_private->m_currentframe;
}

void Molecule::SetCurrentFrame ( size_t i )
{
    m_private->m_currentframe=i;
}

const std::vector<Bond>& Molecule::Bonds() const
{
    return m_private->m_bonds;
}

std::vector<Bond>& Molecule::Bonds()
{
    return m_private->m_bonds;
}

size_t Molecule::NBonds ( size_t i, size_t j ) const
{
    if ( i == j ) return 0;
    int nbonds=0;

    std::vector<size_t> searched;
    searched.push_back ( i );
    bool b=false;
    for ( ;; )
    {
        ++nbonds;
        b=FindInShell ( searched,j );
        if ( b ) break;
        if ( searched.size() == Atoms().size() ) break;

    }

    if ( b ) return nbonds;
    else return 0;

}

bool Molecule::FindInShell ( std::vector<size_t>& searched,size_t j ) const
{

    std::vector<size_t> v=searched;

    for ( std::vector<size_t>::const_iterator at=searched.begin();at!=searched.end();++at )
    {
        std::vector<size_t> nh=Neighbours ( *at );
        v.insert ( v.end(),nh.begin(),nh.end() );
    }
    for ( std::vector<size_t>::const_iterator at=v.begin();at!=v.end();++at )
    {
        if ( std::find ( searched.begin(),searched.end(),*at ) == searched.end() )
            searched.push_back ( *at );
    }
    return ( std::find ( searched.begin(),searched.end(),j ) != searched.end() );
}

const std::vector<PDBResidue*>& Molecule::Residues() const
{
    return m_private->m_residues;
}

std::vector<PDBResidue*>& Molecule::Residues()
{
    return m_private->m_residues;
}


std::list<int> Molecule::Elements() const
{
    std::list<int> nuclei;
    std::vector<Atom>::const_iterator it;
    for ( it=Atoms().begin();it!=Atoms().end();++it )
    {
        if ( std::find ( nuclei.begin(),nuclei.end(),it->Z() ) == nuclei.end() )
            nuclei.push_back ( it->Z() );


    }
    nuclei.sort();
    return nuclei;
}

std::list<std::string> Molecule::ElementSymbols() const
{
    std::list<int> el=Elements();
    std::list<std::string> symbols;
    const PeriodicTable& table=GetPeriodicTable();
    std::map<const std::string,const Atom::atomdata>::const_iterator it;

    for ( std::list<int>::const_iterator at=el.begin();at!=el.end();++at )
    {
        for ( it=table.begin();it!=table.end();++it )
        {
            if ( it->second._Z == ( *at ) )
                symbols.push_back ( it->first );
        }
    }

    return symbols;
}

size_t Molecule::IndexFromPDB(const std::string& pdbname, const std::string& resname, const std::string& resindex) const
{
    size_t index=0;
    for( ;index<Atoms().size();++index)
    {
        const PDBResidue* p=Atoms().at(index).Residue();
        if ( p )
        {
            if ( ( kryomol::toupper(resname) == kryomol::toupper(p->Name())  ) && ( kryomol::toupper(resindex) == kryomol::toupper(p->Index()) ) )
            {
                if ( kryomol::toupper(pdbname ) == kryomol::toupper(Atoms().at(index).PDBName() ) )
                    return index;
            }
        }
    }
    throw kryomol::Exception("Could not get atom index from pdb name");
}


const std::vector<double>& Molecule::Populations() const
{
    return m_private->m_populations;
}

std::vector<double>& Molecule::Populations()
{
    return m_private->m_populations;
}

void Molecule::SuperImpose(size_t refframe)
{
    std::vector<size_t> atoms(Atoms().size());
    for(size_t i=0;i<atoms.size();++i)
        atoms[i]=i;
    SuperImpose(refframe,atoms);
}


void Molecule::SuperImpose(size_t refframe, const std::vector<size_t>& atoms)
{
    if ( atoms.empty() ) throw kryomol::Exception("The atom list is empty");


    for(size_t frame=0;frame<Frames().size();++frame)
    {
        Coordinate centroid(0,0,0);
        
        for(size_t i=0;i<atoms.size();++i)
        {
            centroid.x()+=Frames().at(frame).XYZ().at(atoms[i]).x();
            centroid.y()+=Frames().at(frame).XYZ().at(atoms[i]).y();
            centroid.z()+=Frames().at(frame).XYZ().at(atoms[i]).z();

        }
        centroid.x()/=atoms.size();
        centroid.y()/=atoms.size();
        centroid.z()/=atoms.size();
        
        for(std::vector<Coordinate>::iterator ct=Frames().at(frame).XYZ().begin();ct!=Frames().at(frame).XYZ().end();++ct)
        {
            (*ct)-=centroid;
        }
        
    }

    size_t nrows=3*atoms.size();
    D2Array<double> modelmatrix(nrows,3);
    D1Array<double> solution(3);
    D1Array<double> diff(nrows);
    for(size_t frame=0;frame<Frames().size();++frame)
    {
        size_t nrow=0;
        if ( frame != refframe )
        {
            for(size_t natom=0;natom<atoms.size();++natom)
            {
                size_t atomindex=atoms[natom];
                double resx=Frames().at(frame).XYZ().at(atomindex).x()-Frames().at(refframe).XYZ().at(atomindex).x();
                double resy=Frames().at(frame).XYZ().at(atomindex).y()-Frames().at(refframe).XYZ().at(atomindex).y();
                double resz=Frames().at(frame).XYZ().at(atomindex).z()-Frames().at(refframe).XYZ().at(atomindex).z();

                double xx=       Frames().at(refframe).XYZ().at(atomindex).x()+Frames().at(frame).XYZ().at(atomindex).x();
                double yy= Frames().at(refframe).XYZ().at(atomindex).y()+Frames().at(frame).XYZ().at(atomindex).y();
                double zz=Frames().at(refframe).XYZ().at(atomindex).z()+Frames().at(frame).XYZ().at(atomindex).z();

                modelmatrix(nrow,0)=0;
                modelmatrix(nrow,1)=zz;
                modelmatrix(nrow,2)=-yy;
                diff(nrow)=resx;
                nrow++;
                modelmatrix(nrow,0)=-zz;
                modelmatrix(nrow,1)=0;
                modelmatrix(nrow,2)=xx;
                diff(nrow)=resy;
                nrow++;

                modelmatrix(nrow,0)=yy;
                modelmatrix(nrow,1)=-xx;
                modelmatrix(nrow,2)=0;
                diff(nrow)=resz;
                nrow++;
            }
            //qryomol::SVDFit(diff,modelmatrix,solution);
            for(size_t i=0;i<3;++i)
                std::cout << "sol(i)=" << solution(i) << std::endl;
            double solsum=solution(0)*solution(0)+solution(1)*solution(1)+solution(2)*solution(2);
            //\tan(theta/2)
            double t=sqrt(solsum);
            //obtaint the directo cosines l m and n
            solution(0)/=t;
            solution(1)/=t;
            solution(2)/=t;
            double angle=atan(t);
            std::cout << "angle is " << (180*angle/M_PI) << std::endl;

            double cangle=cos(angle);
            double sangle=sin(angle);

            Quaternion q(cangle,Coordinate(solution(0)*sangle,solution(1)*sangle,solution(2)*sangle));
            Quaternion q1=q.Invert();
            //and now apply rotation
            for(std::vector<Coordinate>::iterator ct=Frames().at(frame).XYZ().begin();ct!=Frames().at(frame).XYZ().end();++ct)
            {
                (*ct)=(q1^(*ct)^q).V();
            }
            std::cout << Frames().at(frame) << std::endl;

        }
    }

}

void Molecule::EckartTransform(size_t refframe, const std::vector<size_t>& atoms)
{
    if ( atoms.empty() ) throw kryomol::Exception("The atom list is empty");

    for(size_t frame=0;frame<Frames().size();++frame)
    {
        Coordinate masscenter(0,0,0);
        double totalmass=0;
        for(size_t i=0;i<atoms.size();++i)
        {
            double mass =this->Atoms().at(i).AtomicMass();
            masscenter.x()+=Frames().at(frame).XYZ().at(atoms[i]).x()*mass;
            masscenter.y()+=Frames().at(frame).XYZ().at(atoms[i]).y()*mass;
            masscenter.z()+=Frames().at(frame).XYZ().at(atoms[i]).z()*mass;
            totalmass+=mass;

        }
        masscenter.x()/=totalmass;
        masscenter.y()/=totalmass;
        masscenter.z()/=totalmass;

        for(std::vector<Coordinate>::iterator ct=Frames().at(frame).XYZ().begin();ct!=Frames().at(frame).XYZ().end();++ct)
        {
            (*ct)-=masscenter;
        }

    }



    size_t nrows=3*atoms.size();
    D2Array<double> modelmatrix(nrows,3);
    D1Array<double> solution(3);
    D1Array<double> diff(nrows);
    for(size_t frame=0;frame<Frames().size();++frame)
    {
        size_t nrow=0;
        if ( frame != refframe )
        {
            for(size_t natom=0;natom<atoms.size();++natom)
            {
                size_t atomindex=atoms[natom];
                double mass =this->Atoms().at(natom).AtomicMass();
                double resx=Frames().at(frame).XYZ().at(atomindex).x()-Frames().at(refframe).XYZ().at(atomindex).x();
                double resy=Frames().at(frame).XYZ().at(atomindex).y()-Frames().at(refframe).XYZ().at(atomindex).y();
                double resz=Frames().at(frame).XYZ().at(atomindex).z()-Frames().at(refframe).XYZ().at(atomindex).z();

                double xx=       Frames().at(refframe).XYZ().at(atomindex).x()+Frames().at(frame).XYZ().at(atomindex).x();
                double yy= Frames().at(refframe).XYZ().at(atomindex).y()+Frames().at(frame).XYZ().at(atomindex).y();
                double zz=Frames().at(refframe).XYZ().at(atomindex).z()+Frames().at(frame).XYZ().at(atomindex).z();

                modelmatrix(nrow,0)=0;
                modelmatrix(nrow,1)=mass*zz;
                modelmatrix(nrow,2)=-mass*yy;
                diff(nrow)=mass*resx;
                nrow++;
                modelmatrix(nrow,0)=-mass*zz;
                modelmatrix(nrow,1)=0;
                modelmatrix(nrow,2)=mass*xx;
                diff(nrow)=mass*resy;
                nrow++;

                modelmatrix(nrow,0)=mass*yy;
                modelmatrix(nrow,1)=-mass*xx;
                modelmatrix(nrow,2)=0;
                diff(nrow)=mass*resz;
                nrow++;
            }
            //qryomol::SVDFit(diff,modelmatrix,solution);
            for(size_t i=0;i<3;++i)
                std::cout << "sol(i)=" << solution(i) << std::endl;
            double solsum=solution(0)*solution(0)+solution(1)*solution(1)+solution(2)*solution(2);
            //\tan(theta/2)
            double t=sqrt(solsum);
            //obtaint the directo cosines l m and n
            solution(0)/=t;
            solution(1)/=t;
            solution(2)/=t;
            double angle=atan(t);
            std::cout << "angle is " << (180*angle/M_PI) << std::endl;

            double cangle=cos(angle);
            double sangle=sin(angle);

            Quaternion q(cangle,Coordinate(solution(0)*sangle,solution(1)*sangle,solution(2)*sangle));
            Quaternion q1=q.Invert();
            //and now apply rotation
            for(std::vector<Coordinate>::iterator ct=Frames().at(frame).XYZ().begin();ct!=Frames().at(frame).XYZ().end();++ct)
            {
                (*ct)=(q1^(*ct)^q).V();
            }
            std::cout << Frames().at(frame) << std::endl;

        }
    }

}


void Molecule::CalculateMassCenter(bool real)
{

    for(size_t frame=0;frame<Frames().size();++frame)
    {
        Coordinate cm(0,0,0);

        for(size_t i=0;i<Atoms().size();++i)
        {
            cm.x()+=Frames().at(frame).XYZ().at(i).x();
            cm.y()+=Frames().at(frame).XYZ().at(i).y();
            cm.z()+=Frames().at(frame).XYZ().at(i).z();

        }
        cm.x()/=Atoms().size();
        cm.y()/=Atoms().size();
        cm.z()/=Atoms().size();

        for(std::vector<Coordinate>::iterator ct=Frames().at(frame).XYZ().begin();ct!=Frames().at(frame).XYZ().end();++ct)
        {
            (*ct)-=cm;
        }

    }
}


/*void Molecule::SetDihedral(size_t i, size_t j, size_t k, size_t l,float dihedral)
{
    float oldd = Coordinate::Dihedral(CurrentFrame().XYZ().at(i),CurrentFrame().XYZ().at(j),CurrentFrame().XYZ().at(k),CurrentFrame().XYZ().at(l));
    float pass = dihedral*M_PI/180.-oldd;

    //m_rotaxis[0]=CurrentFrame().XYZ().at(j);
    //m_rotaxis[1]=CurrentFrame().XYZ().at(k);
    //RotateNeighbours(j,k,1,pass);

    //m_rotatedatoms.clear();

}*/

const float rotbondpass=0.3f;
void Molecule::RotateBond(size_t i, size_t j, int sense)
{
    bool clockwise = ( sense <= 0);
    qDebug() << "sense=" << clockwise;
    this->CurrentFrame().RotateBond(i,j,rotbondpass, clockwise );
    /*m_rotaxis[0]=CurrentFrame().XYZ().at(i);
  m_rotaxis[1]=CurrentFrame().XYZ().at(j);
  float pass=rotbondpass;
  if (sense < 0 ) pass*=-1;
  RotateNeighbours(i,j,sense,pass);
  m_rotatedatoms.clear();*/

}

void Molecule::RotateNeighbours(size_t i, size_t j, int sense,double pass)
{
    //rotate atom
    //std::cout << "Rotate[i,j]=" << i  << " " << j << std::endl;
    qDebug() << m_rotaxis[0].x() << "," << m_rotaxis[1].y();
    std::vector<size_t>::iterator pos=std::find( m_rotatedatoms.begin(), m_rotatedatoms.end() , j );
    if( pos == m_rotatedatoms.end() )
    {
        Coordinate::RotAroundAxis(CurrentFrame().XYZ()[j],m_rotaxis[0],m_rotaxis[1],pass);
        m_rotatedatoms.push_back(j);
    }
    else return;

    std::vector<Bond>::iterator bt;
    for ( bt=Bonds().begin();bt!=Bonds().end();++bt )
    {
        //std::cout << bt->I() << "," << bt->J() << std::endl;
        if ( ( bt->I() == j ) && ( bt->J() != i)   )
        {
            qDebug() << "new bond="<< bt->I() << "," << bt->J();
            qDebug() << "pass" << pass;
            RotateNeighbours( bt->J(), j,sense,pass);
        }

        if (  ( bt->J() == j ) && ( bt->I() != i)   )
        {
            qDebug() << "new bond="<< bt->I() << "," << bt->J();
            qDebug() << "pass" << pass;
            RotateNeighbours( bt->I(), j,sense,pass);
        }

    }
}


float Molecule::GetAnglePlanes(std::vector<size_t>& v, bool degrees)
{
    Coordinate v1=CurrentFrame().XYZ().at(v[2])-CurrentFrame().XYZ().at(v[1]);
    Coordinate v2=CurrentFrame().XYZ().at(v[1])-CurrentFrame().XYZ().at(v[0]);

    Coordinate v3=v1^v2;

    v1=CurrentFrame().XYZ().at(v[5])-CurrentFrame().XYZ().at(v[4]);
    v2=CurrentFrame().XYZ().at(v[4])-CurrentFrame().XYZ().at(v[3]);

    Coordinate v4=v1^v2;

    //Ok lets get the angle between vectors
    float angle=Coordinate::Angle(v3,v4);

    if(degrees)
        return angle*180/M_PI;
    else
        return angle;
}

void Molecule::ResetSelection()
{

    std::vector<Atom>::iterator it;

    for(it=Atoms().begin();it!=Atoms().end();it++)
    {
        it->Select(false);
    }
}

unsigned int Molecule::CountHeavyAtoms() const
{
    unsigned int counter=0;
    for(size_t at=0;at<Atoms().size();++at)
    {
        if ( Atoms().at(at).Symbol() != "H" )  counter++;
    }
    return counter;
}

std::string Molecule::GetEnergyLevel() const
{
    return m_private->m_energylevel;
}

void Molecule::SetEnergyLevel(std::string level) const
{
    m_private->m_energylevel=level;
}


