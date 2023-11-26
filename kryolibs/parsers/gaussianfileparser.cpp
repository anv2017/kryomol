/*****************************************************************************************
                            gaussianfileparser.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "gaussianfileparser.h"
#include "stringtools.h"
#include "molecule.h"
#include "exception.h"
#include "stdlib.h"
#include "mathtools.h"
#include "orbitalarray.h"

#include <QProgressDialog>
#include <QDebug>

#include <vector>
#include <iostream>
#include <cmath>
#include <clocale>
#include <algorithm>

using namespace kryomol;

GaussianFileParser::GaussianFileParser ( const char* file ) : Parser ( file )
{}

GaussianFileParser::GaussianFileParser ( std::istream* stream ) : Parser ( stream )
{}

GaussianFileParser::~GaussianFileParser()
{}

bool GaussianFileParser::ParseFile ( std::streampos pos )
{
    //clear the position vector
    m_pos.clear();

    m_file->clear();
    m_file->seekg ( pos,std::ios::beg );

    m_version=GetVersion();

    m_file->clear();
    m_file->seekg ( pos,std::ios::beg );

    m_level=GetLevel();

    m_file->clear();
    m_file->seekg ( pos,std::ios::beg );

    GetGeometry();

    m_file->clear();
    m_file->seekg ( pos,std::ios::beg );

    if (ExistOrbitals())
    {
        ParseOrbitals (pos);
    }

    m_file->clear();
    m_file->seekg ( pos,std::ios::beg );

    return true;

}

bool GaussianFileParser::GetGeometry()
{
    std::string line;
    bool bfound=false;
    //we should have now the header skip one line
    std::getline ( *m_file,line );
    std::streampos initpos=m_file->tellg();
    bool bscan=false;
    while ( std::getline ( *m_file,line ) )
    {
        if ( !bscan )
            if ( ( line.find ( "Scan                            !" ) != std::string::npos ) //gaussian03
                 || ( line.find ( "Scan                         !" ) != std::string::npos ) ) bscan=true; //gaussian98

        if ( line.find ( "Stationary point found" ) != std::string::npos && !bscan ) break;
        StringTokenizer token ( line );
        if ( token.size() )
        {
            if ( line.find ( "Standard orientation:" ) != std::string::npos )
            {
                bfound=true;

                for ( int i=0;i<3;i++ )
                {
                    std::getline ( *m_file,line );

                }

                m_pos.push_back ( m_file->tellg() );
            }

        }

    }
    //Kein Standard Orientation, try input
    if ( m_pos.empty() )
    {
        m_file->clear();
        m_file->seekg ( initpos,std::ios::beg );
        while ( std::getline ( *m_file,line ) )
        {
            if ( line.find ( "Input orientation:" ) != std::string::npos )
            {
                bfound=true;

                for ( int i=0;i<3;i++ )
                {
                    std::getline ( *m_file,line );
                }

                m_pos.push_back ( m_file->tellg() );
            }
            if ( !line.empty() )
                if ( line.at(0)  == '#'  ) break; //split different jobs
        }

    }
    //Kein Standard Orientation or Input Orientation, try Z-matrix
    if ( m_pos.empty() )
    {
        m_file->clear();
        m_file->seekg ( initpos,std::ios::beg );
        while ( std::getline ( *m_file,line ) )
        {
            if ( line.find ( "Z-Matrix orientation:" ) != std::string::npos )
            {
                bfound=true;

                for ( int i=0;i<3;i++ )
                {
                    std::getline ( *m_file,line );
                }

                m_pos.push_back ( m_file->tellg() );
            }
            if ( !line.empty() )
                if ( line.at(0)  == '#'  ) break; //split different jobs
        }

    }
    if ( m_pos.empty() ) return false;

    std::vector<std::streampos>::iterator pt;
    bool bthreshold=false;
    Threshold threshold;
    //Get simply the last molecule
    Molecules()->push_back ( Molecule() );
    pt=m_pos.begin();
    m_file->clear();
    m_file->seekg ( *pt,std::ios::beg );
    std::getline ( *m_file,line ); //skip one line
    while ( std::getline ( *m_file,line ) )
    {
        StringTokenizer token ( line," " );
        size_t size=token.size();

        if ( size >= 5 )
        {
            Molecules()->back().Atoms().push_back ( Atom ( atoi ( token.at ( 1 ).c_str() ) ) );
        }
        else break;
    }
    for ( pt=m_pos.begin();pt!=m_pos.end();pt++ )
    {
        Molecules()->back().Frames().push_back ( Frame ( &Molecules()->back() ) );

        m_file->clear();
        m_file->seekg ( *pt,std::ios::beg );
        std::getline ( *m_file,line ); //skip one line
        while ( std::getline ( *m_file,line ) )
        {
            StringTokenizer token ( line," " );
            size_t size=token.size();

            if ( size >= 5 )
            {
                Coordinate c;
                c.x() =std::atof ( token.at ( size-3 ).c_str() );
                c.y() =std::atof ( token.at ( size-2 ).c_str() );
                c.z() =std::atof ( token.at ( size-1 ).c_str() );
                Molecules()->back().Frames().back().XYZ().push_back ( c );
            }

            else  break;
        }
        bool bbreak=false;
        while ( std::getline ( *m_file,line ) )
        {
            if ( line.find ( "Stationary point found" ) != std::string::npos ) break;
            StringTokenizer token ( line," \t," );
            if ( token.size() >= 4 )
            {
                switch ( m_level )
                {
                case MNDO:
                    std::cout << "MNDO" << std::endl;
                    if ( m_version != GAUSSIAN09 )
                    {
                        if ( line.find ( "Energy=" ) != std::string::npos )
                        {
                            Molecules()->back().Frames().back().SetEnergy ( atof ( token.at ( 1 ) ),"SCF" );
                        }
                    }
                    else
                    {
                        if ( token.size() >= 5 )
                            if ( token.at ( 0 ) =="SCF" && token.at ( 1 ) =="Done:" )
                            {
                                Molecules()->back().Frames().back().SetEnergy ( std::atof ( token.at ( 4 ).c_str() ),"SCF" );
                            }

                    }
                    break;
                case MP2:
                    if ( token.size() >= 5 )
                        if ( token.at ( 3 ).find ( "MP2" ) != std::string::npos )
                            Molecules()->back().Frames().back().SetEnergy ( atof ( token.at ( 5 ) ),"MP2" );
                    break;
                case QCISD:

                    if ( token[2]=="E(CORR)=" )
                    {
                        Molecules()->back().Frames().back().SetEnergy ( atof ( token.at ( 3 ) ),"QCISD" );
                    }
                    break;
                case SCRF:
                    if ( m_version != GAUSSIAN09 )
                    {
                        if ( line.find ( "with all non electrostatic terms" ) != std::string::npos )
                            Molecules()->back().Frames().back().SetEnergy ( atof ( token.back() ),"SCRF" );
                        if ( line.find ( "Total energy (include solvent energy)" ) != std::string::npos )
                            Molecules()->back().Frames().back().SetEnergy ( atof ( token.back() ),"SCRF" );
                    }
                    else
                    {
                        if ( token.size() >= 5 )
                            if ( token.at ( 0 ) =="SCF" && token.at ( 1 ) =="Done:" )
                            {
                                Molecules()->back().Frames().back().SetEnergy ( std::atof ( token.at ( 4 ).c_str() ),"SCF" );

                            }
                    }
                    break;
                case HF:
                default:
                    if ( token.size() >= 5 )
                        if ( token.at ( 0 ) =="SCF" && token.at ( 1 ) =="Done:" )
                        {
                            Molecules()->back().Frames().back().SetEnergy ( std::atof ( token.at ( 4 ).c_str() ),"SCF" );

                        }
                    break;
                }

                if ( line.find ( "before annihilation" ) !=std::string::npos )
                {
                    Molecules()->back().Frames().back().SetS2 ( std::atof ( token.at ( 3 ).c_str() ) );

                }
                if ( line.find ( "Forces (Hartrees/Bohr)" ) !=std::string::npos ) GetForces ( Molecules()->back() );

                if ( token.at ( 0 ) =="Maximum" && token.at ( 1 ) =="Force" )
                {
                    Molecules()->back().Frames().back().SetMaximumForce ( std::atof ( token.at ( 2 ).c_str() ) );
                    if ( !bthreshold )
                        threshold.maxforce = std::atof ( token.at ( 3 ).c_str() );

                }

                if ( token.at ( 0 ) =="RMS" && token.at ( 1 ) =="Force" )
                {
                    Molecules()->back().Frames().back().SetRMSForce ( std::atof ( token.at ( 2 ).c_str() ) );
                    std::cout << "RMSForce: " << std::atof ( token.at ( 2 ).c_str() ) << " " << Molecules()->back().Frames().back().GetRMSForce() << std::endl;
                    if ( !bthreshold )
                    {
                        threshold.rmsforce = std::atof ( token.at ( 3 ).c_str() );
                        std::cout << "Threshold RMSForce: " << std::atof ( token.at ( 3 ).c_str() ) << std::endl;
                    }
                }

                if ( token.at ( 0 ) =="Maximum" && token.at ( 1 ) =="Displacement" )
                {
                    Molecules()->back().Frames().back().SetMaximumDisplacement ( std::atof ( token.at ( 2 ).c_str() ) );
                    if ( !bthreshold )
                        threshold.maxdisplacement = std::atof ( token.at ( 3 ).c_str() );

                }


                if ( token.at ( 0 ) =="RMS" && token.at ( 1 ) =="Displacement" )
                {
                    Molecules()->back().Frames().back().SetRMSDisplacement ( std::atof ( token.at ( 2 ).c_str() ) );
                    if ( !bthreshold )
                        threshold.rmsdisplacement = std::atof ( token.at ( 3 ).c_str() );
                    bthreshold=true;
                    bbreak=true;

                }


            }
            if ( bbreak ) break;

        }


        Molecules()->back().Frames().back().SetThreshold ( threshold );

        if ( pt+1 != m_pos.end() )
        {
            GetDipole ( *pt,* ( pt+1 ) );
            GetESPCharges ( *pt,* ( pt+1 ) );
        }
        else
        {
            m_file->seekg ( 0,std::ios::end );
            std::streampos end=m_file->tellg();
            GetDipole ( *pt,end );
            GetESPCharges ( *pt,end );
        }
    }

    Molecules()->back().SetBonds();

    return true;
}


bool GaussianFileParser::ParseOrbitals(std::streampos pos)
{
    m_file->clear();
    m_file->seekg ( pos,std::ios::beg );
    bool b = false;

    m_beta = ExistAlphaBetaOrbitals();

    m_file->clear();
    m_file->seekg ( pos,std::ios::beg );
    if (GetCoordinatesType())
    {
        m_file->clear();
        m_file->seekg ( pos,std::ios::beg );
        if (GetHomoLumo())
        {
            m_file->clear();
            m_file->seekg ( pos,std::ios::beg );
            if (GetBasisCenters())
            {
                m_file->clear();
                m_file->seekg ( pos,std::ios::beg );
                if (m_beta)
                    b = GetAlphaBetaOrbitalData();
                else
                    b = GetOrbitalData();
            }
        }
    }
    return b;
}

bool GaussianFileParser::GetBasisCenters()
{
    std::string line;
    while ( std::getline ( *m_file,line ) )
    {
        if ( line.find ( "AO basis set in the form of general basis input" ) != std::string::npos )
        {
            m_norbitals = 0;
            std::getline( *m_file,line );
            StringTokenizer token ( line," " );
            size_t sizeAtom=token.size();
            while (sizeAtom == 2)
            {
                std::vector<Orbital> orbitals;

                StringTokenizer tokenAtom ( line," " );
                int atom = atof (tokenAtom.at(0))-1;

                std::getline( *m_file,line );
                StringTokenizer token ( line," \t\r" );
                size_t sizeOrbital = token.size();
                while (sizeOrbital==4)
                {
                    StringTokenizer tokenOrbital ( line," \t\r" );
                    std::vector<float> xs;
                    std::vector<float> xp;
                    std::vector<float> alpha;

                    for (int i=0; i<atof(tokenOrbital.at(1)); i++)
                    {

                        std::getline( *m_file,line );
                        StringTokenizer tokenCoefficient ( line," \t\r" );
                        if (tokenOrbital.at(0)=="SP")
                        {
                            if (tokenCoefficient.size()>=3)
                            {
                                xs.push_back(atof(tokenCoefficient.at(1)));
                                xp.push_back(atof(tokenCoefficient.at(2)));
                                alpha.push_back(atof(tokenCoefficient.at(0)));
                            }
                        }
                        else
                        {
                            if (tokenCoefficient.size()>=2)
                            {
                                xs.push_back(atof(tokenCoefficient.at(1)));
                                alpha.push_back(atof(tokenCoefficient.at(0)));
                            }
                        }
                    }

                    if (tokenOrbital.at(0)=="S")
                    {
                        orbitals.push_back( Orbital (Orbital::S,alpha,xs,xp));
                        m_norbitals += 1;
                    }
                    if (tokenOrbital.at(0)=="SP")
                    {
                        orbitals.push_back( Orbital (Orbital::SP,alpha,xs,xp));
                        m_norbitals += 4;
                    }
                    if (tokenOrbital.at(0)=="P")
                    {
                        orbitals.push_back(Orbital (Orbital::P,alpha,xs,xp));
                        m_norbitals += 3;
                    }
                    if (tokenOrbital.at(0)=="D")
                    {
                        orbitals.push_back(Orbital (Orbital::D,alpha,xs,xp));
                        m_norbitals += m_typeD;
                    }
                    if (tokenOrbital.at(0)=="F")
                    {
                        orbitals.push_back(Orbital (Orbital::F,alpha,xs,xp));
                        m_norbitals += m_typeF;
                    }

                    std::getline( *m_file,line );
                    StringTokenizer token ( line," \t\r" );
                    sizeOrbital = token.size();
                }

                std::getline( *m_file,line );
                StringTokenizer tokenAux ( line," " );
                sizeAtom=tokenAux.size();

                m_atoms.push_back(atom);
                m_orbitals.push_back(orbitals);
            }

            return true;
        }
    }
    return false;
}

bool GaussianFileParser::GetOrbitalData()
{
    bool b =false;
    int frame =0;
    std::vector<Frame>::iterator ft = Molecules()->back().Frames().begin();
    size_t pt=0;
    std::string line;
    while ( std::getline(*m_file,line) )
    {
        if ( line.find ( "Molecular Orbital Coefficients" ) != std::string::npos )
        {
            while (!(pt==m_pos.size()-1))
            {
                if ( (m_file->tellg()>m_pos.at(pt)) && (m_file->tellg()<m_pos.at(pt+1)) )
                    break;
                else
                {
                    ++pt;
                    ++ft;
                    frame+=1;
                }
            }

            D2Array<float> matrix;
            std::vector<float> eigenvalues;
            std::vector<float> occupations;
            matrix.Initialize(m_norbitals,m_norbitals,0.0);

            std::getline( *m_file,line );
            size_t K=0;
            while ( (line.find("Density Matrix:") == std::string::npos ) && (line.find("Condensed to atoms (all electrons):") == std::string::npos ) )
            {
                StringTokenizer tokenSize ( line," " );
                size_t M=tokenSize.size();

                std::getline( *m_file,line );
                //get the occupation
                StringTokenizer ocstr(line," ");
                for(StringTokenizer::const_iterator it=ocstr.begin();it!=ocstr.end();++it)
                {
                    if ( it->back() == 'O' )
                    {
                        occupations.push_back(2.0);
                    }
                    else {
                        occupations.push_back(0.0);
                    }

                }
                std::getline( *m_file,line );
                StringTokenizer tokenEigenv ( line," " );
                for (size_t e=0; e<M; e++)
                {
                    eigenvalues.push_back(atof(tokenEigenv.at(tokenEigenv.size()-M+e)));
                }
                for (size_t n=0; n<(size_t)m_norbitals; n++)
                {
                    std::getline( *m_file,line );
                    StringTokenizer tokenCoeff ( line," " );
                    for (size_t m=0; m<M; m++)
                    {
                        matrix(n,K+m) = atof(tokenCoeff.at(tokenCoeff.size()-M+m));
                    }
                }
                K += M;
                std::getline( *m_file,line );
            }
            std::vector<BasisCenter> basis;
            for (size_t i=0; i<m_atoms.size(); ++i)
                basis.push_back(BasisCenter(ft->XYZ().at(m_atoms.at(i)),m_orbitals.at(i)));
            ft->OrbitalsData().SetHomo(m_homo);
            ft->OrbitalsData().SetLumo(m_lumo);
            ft->OrbitalsData().SetTypeD(m_typeD);
            ft->OrbitalsData().SetTypeF(m_typeF);
            ft->OrbitalsData().SetBasisCenters(basis);
            ft->OrbitalsData().SetCoefficients(matrix);
            ft->OrbitalsData().SetEigenvalues(eigenvalues);
            ft->OrbitalsData().SetOccupations(occupations);
            b = true;

            ++pt;
            ++ft;
            frame+=1;
        }
    }
    return b;
}

bool GaussianFileParser::GetAlphaBetaOrbitalData()
{
    bool b =false;
    int frame =0;
    std::vector<Frame>::iterator ft = Molecules()->back().Frames().begin();
    size_t pt=0;
    std::string line;
    while ( std::getline(*m_file,line) )
    {
        if ( line.find ( "Molecular Orbital Coefficients" ) != std::string::npos )
        {
            while (!(pt==m_pos.size()-1))
            {
                if ( (m_file->tellg()>m_pos.at(pt)) && (m_file->tellg()<m_pos.at(pt+1)) )
                    break;
                else
                {
                    ++pt;
                    ++ft;
                    frame+=1;
                }
            }

            D2Array<float> matrix;
            D2Array<float> betamatrix;
            std::vector<float> eigenvalues;
            std::vector<float> betaeigenvalues;
            matrix.Initialize(m_norbitals,m_norbitals,0.0);
            betamatrix.Initialize(m_norbitals,m_norbitals,0.0);

            if ( line.find ( "Alpha Molecular Orbital Coefficients:" ) != std::string::npos )
            {
                size_t K = 0;

                std::getline( *m_file,line );
                while((line.find("Beta Molecular Orbital Coefficients:")) == std::string::npos )
                {
                    StringTokenizer tokenSize ( line," " );
                    size_t M=tokenSize.size();

                    std::getline( *m_file,line );
                    std::getline( *m_file,line );
                    StringTokenizer tokenEigenv ( line," " );
                    for (size_t e=0; e<M; e++)
                    {
                        eigenvalues.push_back(atof(tokenEigenv.at(tokenEigenv.size()-M+e)));
                    }
                    for (size_t n=0; n<(size_t)m_norbitals; n++)
                    {
                        std::getline( *m_file,line );
                        StringTokenizer tokenCoeff ( line," " );
                        for (size_t m=0; m<M; m++)
                        {
                            matrix(n,K+m) = atof(tokenCoeff.at(tokenCoeff.size()-M+m));
                        }
                    }
                    K += M;
                    std::getline( *m_file,line );
                }

                K = 0;
                std::getline( *m_file,line );
                while((line.find("Alpha Density Matrix:")) == std::string::npos )
                {
                    StringTokenizer tokenSize ( line," " );
                    size_t M=tokenSize.size();

                    std::getline( *m_file,line );
                    std::getline( *m_file,line );
                    StringTokenizer tokenEigenv ( line," " );
                    for (size_t e=0; e<M; e++)
                    {
                        betaeigenvalues.push_back(atof(tokenEigenv.at(tokenEigenv.size()-M+e)));
                    }
                    for (size_t n=0; n<(size_t)m_norbitals; n++)
                    {
                        std::getline( *m_file,line );
                        StringTokenizer tokenCoeff ( line," " );
                        for (size_t m=0; m<M; m++)
                        {
                            betamatrix(n,K+m) = atof(tokenCoeff.at(tokenCoeff.size()-M+m));
                        }
                    }
                    K += M;
                    std::getline( *m_file,line );
                }
            }

            std::vector<BasisCenter> basis;
            for (size_t i=0; i<m_atoms.size(); ++i)
                basis.push_back(BasisCenter(ft->XYZ().at(m_atoms.at(i)),m_orbitals.at(i)));
            ft->OrbitalsData().SetHomo(m_homo);
            ft->OrbitalsData().SetLumo(m_lumo);
            ft->OrbitalsData().SetTypeD(m_typeD);
            ft->OrbitalsData().SetTypeF(m_typeF);
            ft->OrbitalsData().SetBasisCenters(basis);
            ft->OrbitalsData().SetCoefficients(matrix);
            ft->OrbitalsData().SetEigenvalues(eigenvalues);
            ft->OrbitalsData().SetBetaCoefficients(betamatrix);
            ft->OrbitalsData().SetBetaEigenvalues(betaeigenvalues);
            b = true;
            ++pt;
            ++ft;
            frame+=1;
        }
    }
    return b;
}

bool GaussianFileParser::ExistOrbitals()
{
    std::string line;
    while ( std::getline ( *m_file,line ))
    {
        if ( line.find("Molecular Orbital Coefficients") != std::string::npos )
            return true;
    }
    return false;
}

bool GaussianFileParser::ExistAlphaBetaOrbitals()
{
    std::string line;
    bool bfound = false;
    while (( std::getline ( *m_file,line ) && !(bfound)))
    {
        if ( line.find("Alpha Orbitals:") != std::string::npos )
            bfound = true;
    }
    return bfound;
}


bool GaussianFileParser::GetCoordinatesType()
{
    std::string line;
    while (std::getline(*m_file,line))
    {
        if (line.find("Standard basis:") != std::string::npos || ( line.find("General basis read from cards:") != std::string::npos ) )
        {
            if ( line.find ( "(5D, 7F)" ) != std::string::npos )
            {
                m_typeD = 5;
                m_typeF = 7;
            }
            if ( line.find ( "(6D, 7F)" ) != std::string::npos )
            {
                m_typeD = 6;
                m_typeF = 7;
            }
            if ( line.find ( "(5D, 10F)" ) != std::string::npos )
            {
                m_typeD = 5;
                m_typeF = 10;
            }
            if ( line.find ( "(6D, 10F)" ) != std::string::npos )
            {
                m_typeD = 6;
                m_typeF = 10;
            }
            return true;
        }
    }
    return false;
}

bool GaussianFileParser::GetHomoLumo()
{
    std::string line;
    while ( std::getline ( *m_file,line ) )
    {
        if ( line.find ( "alpha electrons" ) != std::string::npos )
        {
            StringTokenizer token(line," \t\r");
            m_homo = atof(token.at(0));
            if (m_beta)
                m_lumo = m_homo;
            else
                m_lumo = m_homo+1;

            return true;
        }

        if ( line.find ( "NAE=" ) != std::string::npos )
        {
            StringTokenizer token(line," \t\r");
            for (size_t i=0; i<token.size()-1; i++)
            {
                if (token.at(i) == "NAE=")
                {
                    m_homo = atof(token.at(i+1));
                    if (m_beta)
                        m_lumo = m_homo;
                    else
                        m_lumo = m_homo+1;

                }
            }
            return true;
        }
    }
    return false;
}


void GaussianFileParser::GetTransitionChanges()
{

    std::string line;
    int transition = 0;
    std::vector< std::vector<TransitionChange> > vectortransitions;
    while ( std::getline ( *m_file,line ) )
    {
        if ( line.find ( "Excited State" ) != std::string::npos )
        {
            transition += 1;

            std::vector<TransitionChange> transitions;

            std::getline ( *m_file,line );
            if (m_beta)
            {
                while ((line.find ( " -> " ) != std::string::npos )||(line.find ( " ->" ) != std::string::npos )||(line.find ( "->" ) != std::string::npos ))
                {
                    if (line.find ( " -> " ) != std::string::npos )
                    {
                        StringTokenizer token(line," \t\r");
                        float coeff = atof(token.at(3));
                        transitions.push_back( TransitionChange(token.at(0),token.at(2),coeff) );
                    }
                    else
                    {
                        if (line.find ( " ->" ) != std::string::npos )
                        {
                            StringTokenizer token1(line,">");
                            StringTokenizer token2(token1.at(0)," ");
                            int i = atof(token2.at(0));
                            StringTokenizer token(token1.at(1)," \t\r");
                            int j = atof(token.at(0));
                            float coeff = atof(token.at(1));
                            transitions.push_back( TransitionChange(i,j,coeff) );
                        }
                        else
                        {
                            StringTokenizer token1(line,">");
                            std::string s = token1.at(0).substr(0,s.size()-2);
                            int i = atof(s);
                            StringTokenizer token(token1.at(1)," \t\r");
                            int j = atof(token.at(0));
                            float coeff = atof(token.at(1));
                            transitions.push_back( TransitionChange(i,j,coeff) );
                        }
                    }

                    std::getline ( *m_file,line );
                }
            }
            else
            {
                while (line.find ( "->" ) != std::string::npos )
                {
                    StringTokenizer tokenstring(line,">");
                    int i = atof(tokenstring.at(0));
                    StringTokenizer token(tokenstring.at(1)," \t\r");
                    int j = atof(token.at(0));
                    float coeff = atof(token.at(1));

                    transitions.push_back( TransitionChange(i,j,coeff) );

                    std::getline ( *m_file,line );
                }
            }
            vectortransitions.push_back(transitions);
        }
    }

    Molecules()->back().Frames().back().SetTransitionChanges(vectortransitions);
}


bool GaussianFileParser::GetESPCharges ( const std::streampos& begin, const std::streampos& end )
{
    m_file->clear ( );
    m_file->seekg ( begin, std::ios::beg );   // moves pointer to

    std::string line;
    std::streampos pos=0;
    while ( std::getline ( *m_file,line ) )
    {
        if ( line.find ( "Charges from ESP fit" ) != std::string::npos )
            break;
    }
    pos=m_file->tellg();

    if ( pos <= 0 ) return false;
    if ( end > 0 && pos > end ) return false;

    Molecules()->back().Frames().back().Charges (Frame::ESP)->reserve ( Molecules()->back().Atoms().size() );
    std::getline ( *m_file,line );
    std::getline ( *m_file,line );
    while ( std::getline ( *m_file,line ) )
    {
        StringTokenizer tok ( line );
        if ( tok.size() >= 3 )
        {
            Molecules()->back().Frames().back().Charges ( Frame::ESP )->push_back ( std::atof ( tok.at ( 2 ).c_str() ) );
        }
        if ( line.find ( "--" ) != std::string::npos ) break;
    }
    return true;
}


void GaussianFileParser::GetForces ( Molecule& molecule )
{
    std::string line;
    molecule.Frames().back().Gradient().reserve(molecule.Atoms().size());
    for ( int i=0;i<2;i++ ) std::getline ( *m_file,line );
    while ( std::getline ( *m_file,line ) )
    {
        if ( line.find ( "---" ) !=std::string::npos )
        {
            return;
        }
        StringTokenizer token ( line," \t\r" );
        Coordinate c;
        c.x()=std::atof ( token.at ( 2 ).c_str() );
        c.y()=std::atof ( token.at ( 3 ).c_str() );
        c.z()=std::atof ( token.at ( 3 ).c_str() );
        molecule.Frames().back().Gradient().push_back(c);
    }


}



GaussianFileParser::gaussversion GaussianFileParser::GetVersion()
{
    std::string line;

    while(std::getline(*m_file,line))
    {
        if ( line.find("Gaussian 09:") != std::string::npos )
            return  GAUSSIAN09;
    }
    return UNDEFINED;

}

QuantumLevel GaussianFileParser::GetLevel()
{
    std::string line, pline;


    while ( line.find ( "-----" ) == std::string::npos )
    {
        line+=pline;
        std::getline ( *m_file,pline );
        line.erase ( 0,1 );
    }
    kryomol::toupper ( line );
    if ( ( line.find ( "MNDO",0 ) != std::string::npos ) ||
         ( line.find ( "AM1",0 ) != std::string::npos ) ||
         ( line.find ( "PM3",0 ) != std::string::npos ) )
    {
        std::cout << "MNDO calculation" << std::endl;
        return MNDO;
    }
    if ( line.find ( "MP2",0 ) !=std::string::npos )
    {
        std::cout << "MP2 calculation" << std::endl;
        return MP2;
    }
    if ( line.find ( "QCISD",0 ) !=std::string::npos )
    {
        std::cout << "QCISD calculation" << std::endl;
        return QCISD;
    }
    if ( line.find ( "SCRF",0 ) !=std::string::npos )
    {
        //this is always present in frequency calculations
        if ( line.find ( "SCRF=CHECK" )  == std::string::npos )
        {
            std::cout << "SCRF calculation" << std::endl;
            return SCRF;
        }
    }

    //Hartree-Fock by default
    std::cout << "HF/DFT calculation" << std::endl;
    return HF;
}



bool GaussianFileParser::HasKeyword ( std::string& line )
{

    StringTokenizer token ( line," " );

    if ( token.size() ==7 )
    {
        if ( token[6]=="analyticall!" )
        {

            return true;
        }
        else
            return false;
    }
    else
        return false;
}


std::vector<JobHeader>& GaussianFileParser::Jobs()
{
    //lets get the lines beginning with a #
    std::string line;
    m_file->clear();
    m_file->seekg(0);
    std::streampos pos=m_file->tellg();
    std::getline ( *m_file,line );
    do
    {
        //it could be that by chance a # line is inside the arquive entry, so do this little hack
        // to avoid this
        if ( line.find ( "1\\1\\",0 ) == 1 || line.find ( "1|1|",0 ) == 1 )
        {
            std::getline ( *m_file,line );
            std::getline ( *m_file,line );
            std::getline ( *m_file,line );
        }
        if ( line.find ( " #",0 ) == 0 ) //#should be the second character and not -#-
        {

            std::string wholeline;
            //take into account that route can be splitted into several lines
            while ( line.find ( "-----" ) == std::string::npos )
            {
                //trim left whitespace
                line.erase( line.begin(), std::find_if( line.begin(), line.end(),
                                                        std::not1( std::ptr_fun( &::isspace ) ) ) );
                wholeline+=line;
                std::getline ( *m_file,line );
                //line.erase ( 0,1 );

            }
            //delete return carriage
            std::string::size_type sit=0;
            while ( sit != std::string::npos )
            {
                sit=wholeline.find_first_of ( "\n\r",sit );
                if ( sit != std::string::npos )
                    wholeline.erase ( sit,1 );
            }

            std::cout << "route=" << wholeline << std::endl;
            JobType tjob=GetJobFromRoute(wholeline);
            m_jobpos.push_back ( JobHeader( tjob,pos ) );
            do
            {
                std::getline(*m_file,line);
            } while ( line.find ( "-----" ) == std::string::npos );
        }
        pos=m_file->tellg();
    }
    while ( std::getline ( *m_file,line ) );

    return m_jobpos;
}

JobType GaussianFileParser::GetJobFromRoute(const std::string& route)
{
    //capitalize
    JobType t=singlepoint;
    std::string croute=kryomol::toupper(route);
    std::cout << "finding job" << std::endl;
    if ( croute.find("OPT") != std::string::npos  )
    {
        t=opt;
        return t;
    }

    if ( croute.find("ADMP") != std::string::npos || croute.find("BOMD") != std::string::npos )
    {
        t=dyn;
        return t;
    }

    if ( croute.find("FREQ") != std::string::npos )
    {
        t=freq;
        return t;
    }

    if ( croute.find("TD") != std::string::npos )
    {
        t=uv;
        return t;
    }

    if ( croute.find("NMR") != std::string::npos )
    {
        t=nmr;
        return t;
    }
    //single point if nothing matches
    return t;
}
bool GaussianFileParser::ParseFrequencies ( std::streampos pos )
{
    m_file->clear();
    m_file->seekg ( pos,std::ios::beg );
    if ( !GetFrequencies() ) return false;
    m_file->clear();
    m_file->seekg ( pos,std::ios::beg );
    //Ok lets get the vectors;
    scounter=0;
    std::string line;
    Molecule& molecule=Molecules()->back();
    molecule.Frames().back().AllocateVectors();
    while ( std::getline ( *m_file,line ) )
    {
        StringTokenizer token ( line," " );
        if ( token.size() >=2 )
            if ( token.at ( 1 ) =="AN" )
            {
                GetVectors();
            }
    }
    if ( !ParseArquive ( pos ) ) std::cerr << "No arquive entry found for frequency computation" << std::endl;

    std::vector<Frequency> v = Molecules()->back().Frames().back().GetFrequencies();
    qDebug() << "LIST FREQUENCIES in ParseFrequencies: " << v.size() << endl;

    return true;
}


bool GaussianFileParser::ParseArquive ( std::streampos pos )
{

    std::string archive;
    std::string line;
    m_file->clear();
    m_file->seekg ( pos,std::ios::beg );
    std::streampos arpos;
    Molecule& molecule=Molecules()->back();
    while ( std::getline ( *m_file,line ) )
    {
        if ( ( line.find ( "1\\1\\" ) != std::string::npos ) || ( line.find ( "1|1|" ) !=  std::string::npos ) ) break;
        arpos=m_file->tellg();
    }
    //windows
    m_file->seekg ( arpos,std::ios::beg );

    while ( std::getline ( *m_file,line ) )
    {
        archive+=line;
    }
    //OK lets clear bakspaces
    std::string::size_type sit=0;
    while ( sit != std::string::npos )
    {
        sit=archive.find_first_of ( " \n\t\r",sit );
        if ( sit != std::string::npos )
            archive.erase ( sit,1 );
    }

    TokByString t ( archive,"\\\\" );

    if ( t.size() < 4 )
        t=TokByString ( archive,"||" );

    std::string& geometry= t.at ( 3 );
    TokByString::iterator tit;

    if ( t.size() < 4 ) return false;

    StringTokenizer g1 ( geometry,"\\|" );

    std::vector<Coordinate>& iorient= molecule.InputOrientation();
    iorient.resize (molecule.Atoms().size());

    size_t i=0;
    StringTokenizer::iterator it=g1.begin();
    it++;
    for ( i=0;it!=g1.end();it++,i++ )
    {
        StringTokenizer g2 ( *it,"," );
        std::cerr << *it << std::endl;
        StringTokenizer::reverse_iterator it2=g2.rbegin();
        iorient[i].z()= std::atof( ( it2++ )->c_str() );
        iorient[i].y()= std::atof ( ( it2++ )->c_str() );
        iorient[i].x()= std::atof ( it2->c_str() );
    }

    //Get the forces
    //OK I have to had three more blocks to get forces and hessian, first
    // I can skip, second is hessian and third are forces
    if ( t.size() < 7 )
    {

        return false;
    }

    molecule.Frames().back().AllocateForces();
    molecule.Frames().back().AllocateHessian();

    //Lets get forces in seventh block
    D1Array<double>& forces=molecule.Frames().back().GetForces();
    g1=StringTokenizer ( t.at ( 6 ),"," );
    for ( it=g1.begin(),i=0;it!=g1.end();it++,i++ )
    {
        forces ( i ) =std::atof ( it->c_str() );
    }


    D2Array<double>& hessian=molecule.Frames().back().GetHessian();

    StringTokenizer g2 ( t.at ( 5 ) );
    it=g2.begin();

    size_t j;
    for ( i=0;i< 3*molecule.Atoms().size();i++ )
    {
        for ( j=0;j <= i ;it++,j++ )
        {
            hessian ( i,j ) =std::atof ( it->c_str() );
            if ( i!=j ) hessian ( j,i ) =hessian ( i,j );
        }
    }

    std::vector<Frequency> v = Molecules()->back().Frames().back().GetFrequencies();
    qDebug() << "LIST FREQUENCIES in ParseArchive: " << v.size() << endl;

    return true;
}


bool GaussianFileParser::ParseUV ( std::streampos pos )
{
    m_file->clear();
    m_file->seekg ( pos,std::ios::beg );
    std::string line;
    Molecule& molecule=Molecules()->back();
    std::vector<Spectralline>& lines=molecule.Frames().back().GetSpectralLines();
    while ( std::getline ( *m_file,line ) )
    {
        StringTokenizer token ( line," =\t" );
        size_t size=token.size();
        if ( size > 2 )
        {
            if ( token.at ( 0 ) =="Excited" && token.at ( 1 ) == "State" )
            {
                Spectralline line;
                line.x=std::atof ( token.at ( 6 ).c_str() );
                line.y0=std::atof ( token.at ( 9 ).c_str() ); //oscillator strength
                lines.push_back ( line );
            }
        }

    }

    //Get Rotatory strenghts in velocity
    m_file->clear();
    m_file->seekg ( pos,std::ios::beg );
    std::vector<Spectralline>::iterator spt=lines.begin();
    std::streampos lastpos=m_file->tellg();
    std::vector<std::streampos> positions;
    bool oldversion = false;
    while(std::getline(*m_file,line) )
    {
        positions.push_back(m_file->tellg() );
        if ( line.find("Total R(velocity) tensor for State=") != std::string::npos )
        {
            oldversion = true;
            std::streampos oldpos=m_file->tellg();
            m_file->clear();
            m_file->seekg(positions.at(positions.size()-3),std::ios::beg);
            std::getline(*m_file,line);
            StringTokenizer tok(line," t\r");
            spt->SetRotatoryStrengthVelocity(std::atof(tok.back().c_str()));
            ++spt;
            positions.clear();
            m_file->seekg(oldpos,std::ios::beg);
        }
    }
    if (!(oldversion))
    {
        m_file->clear();
        m_file->seekg ( pos,std::ios::beg );

        while ( std::getline ( *m_file,line ) )
        {
            StringTokenizer token ( line );
            size_t size=token.size();
            if ( size >= 5 )
            {
                if ( token.at ( 4 ) == "R(velocity)" )
                {

                    for ( std::vector<Spectralline>::iterator it=lines.begin();it!=lines.end();++it )
                    {
                        std::getline ( *m_file,line );
                        StringTokenizer token ( line );
                        it->SetRotatoryStrengthVelocity(std::atof ( token.at ( 4 ).c_str() ) );

                    }
                }
            }
        }
    }
    //Get Now Rotatory Strengths
    m_file->clear();
    m_file->seekg ( pos,std::ios::beg );

    while ( std::getline ( *m_file,line ) )
    {
        StringTokenizer token ( line );
        size_t size=token.size();
        if ( size >= 5 )
        {
            if ( token.at ( 4 ) =="R(length)" )
            {

                for ( std::vector<Spectralline>::iterator it=lines.begin();it!=lines.end();++it )
                {
                    std::getline ( *m_file,line );
                    StringTokenizer token ( line );
                    it->SetRotatoryStrengthLength(std::atof ( token.at ( 4 ).c_str() ) );

                }
            }
        }
    }
    GetTransitionVectors(pos);

    m_file->clear();
    m_file->seekg ( pos,std::ios::beg );

    GetTransitionChanges();

    return true;
}

void GaussianFileParser::GetTransitionVectors(std::streampos pos)
{
    Molecule& molecule=Molecules()->back();
    std::vector<Spectralline>& lines=molecule.Frames().back().GetSpectralLines();
    m_file->clear();
    m_file->seekg(pos,std::ios::beg);
    std::string line;
    while ( std::getline(*m_file,line))
    {
        if  ( line.find("Ground to excited state transition electric dipole") != std::string::npos )
        {
            std::getline(*m_file,line); //skip one line;
            for(std::vector<Spectralline>::iterator it=lines.begin();it!=lines.end();++it)
            {
                std::getline(*m_file,line);
                StringTokenizer tok(line," \t\r");
                it->ElectricDipole().x()=std::atof(tok.at(1).c_str());
                it->ElectricDipole().y()=std::atof(tok.at(2).c_str());
                it->ElectricDipole().z()=std::atof(tok.at(3).c_str());
            }
            std::getline(*m_file,line); //skip one line;
            std::getline(*m_file,line); //skip one line;
            for(std::vector<Spectralline>::iterator it=lines.begin();it!=lines.end();++it)
            {
                std::getline(*m_file,line);
                StringTokenizer tok(line," \t\r");
                it->VelocityDipole().x()=std::atof(tok.at(1).c_str());
                it->VelocityDipole().y()=std::atof(tok.at(2).c_str());
                it->VelocityDipole().z()=std::atof(tok.at(3).c_str());
            }
            std::getline(*m_file,line); //skip one line;
            std::getline(*m_file,line); //skip one line;
            for(std::vector<Spectralline>::iterator it=lines.begin();it!=lines.end();++it)
            {
                std::getline(*m_file,line);
                StringTokenizer tok(line," \t\r");
                it->MagneticDipole().x()=std::atof(tok.at(1).c_str());
                it->MagneticDipole().y()=std::atof(tok.at(2).c_str());
                it->MagneticDipole().z()=std::atof(tok.at(3).c_str());
            }
            return;
        }
    }
    return;
}

bool GaussianFileParser::GetVectors()
{
    int counter=0;
    std::string line;

    //OK put vectors onto last molecule
    Molecule& molecule=Molecules()->back();

    while ( std::getline ( *m_file,line ) )
    {
        StringTokenizer token ( line," " );;

        if ( token.size() <5 )
        {
            scounter+=3;
            return true;
        }
        else
        {
            int nfreq= ( token.size()-2 ) /3;

            std::vector< std::vector<Coordinate> >& coordinates = molecule.GetMode();

            for ( int i=0;i<nfreq;i++ )
            {
                coordinates.at(i+scounter).at(counter).x() = std::atof ( token[2+3*i].c_str() );
                coordinates.at(i+scounter).at(counter).y() = std::atof ( token[2+3*i+1].c_str() );
                coordinates.at(i+scounter).at(counter).z() = std::atof ( token[2+3*i+2].c_str() );

            }
            counter++;
        }
    }

    return true;

}
bool GaussianFileParser::GetFrequencies()
{
    Molecule& molecule=Molecules()->back();

    std::vector<Frequency>& frequencies=molecule.Frames().back().GetFrequencies();

    std::string line;

    size_t counter=0;
    size_t counter1=0;
    size_t counter2=0;
    while ( std::getline ( *m_file,line ) )
    {
        StringTokenizer tok ( line );
        if ( tok.size() )
        {
            int steps=0;
            if ( tok[0]=="Frequencies" )
            {

                for ( size_t j=2;j<tok.size();j++ )
                {
                    frequencies.push_back ( Frequency ( std::atof ( tok[j].c_str() ),0.0f ) );
                    counter++;
                }

            }

            if ( tok[0]=="IR" )
            {
                for ( size_t k=3;k<tok.size();k++ )
                {
                    frequencies[counter1++].y=std::atof ( tok[k].c_str() );
                    steps++;
                }
            }
            if ( tok[0]=="Rot." )
            {
                for ( size_t k=3;k<tok.size();k++ )
                {
                    frequencies[counter2++].z=std::atof ( tok[k].c_str() );
                }

            }

        }
    }

    qDebug() << "LIST FREQUENCIES in GetFrequencies: " << frequencies.size() << endl;

    return true;

}


void GaussianFileParser::ParseChemicalShifts()
{

    const char* locale=std::setlocale(LC_NUMERIC,"C");
    m_file->clear();
    m_file->seekg(0,std::ios::beg);

    Frame& frame=Molecules()->back().CurrentFrame();
    std::string line;
    bool tensorsfound=false;
    while( std::getline(*m_file,line) )
    {
        if ( line.find("Magnetic shielding") != std::string::npos )
            break;
    }
    while ( std::getline ( *m_file,line ) )
    {
        if ( line.find ( "Isotropic =" ) != std::string::npos )
        {
            if ( !tensorsfound) tensorsfound=true;
            frame.CShiftTensors().push_back ( D2Array<double> ( 3,3 ) );
            D2Array<double>& t=frame.CShiftTensors().back();
            for ( int i=0;i<3;++i )
            {
                std::getline ( *m_file,line );
                StringTokenizer token ( line," \t\r" );
                t ( 0,i ) =std::atof ( token.at ( 1 ).c_str() );
                t ( 1,i ) =std::atof ( token.at ( 3 ).c_str() );
                t ( 2,i ) = std::atof ( token.at ( 5 ).c_str() );
            }

        }

    }
    if ( !tensorsfound )
    {
        throw kryomol::Exception("Chemical Shift Tensors not found" );
    }

    std::setlocale(LC_NUMERIC,locale);
}



void GaussianFileParser::ParseCouplingConstants(std::vector<QuantumCoupling>& c)
{
    const char* locale=std::setlocale(LC_NUMERIC,"C");
    std::string line;
    m_file->clear();
    m_file->seekg(0);
    while(std::getline(*m_file,line) )
    {
        if ( line.find("Total nuclear spin-spin coupling J") != std::string::npos )
        {
            while(ExtractJBlock(c))
            {};
            std::setlocale(LC_NUMERIC,locale);
            return;
        }
    }
    std::setlocale(LC_NUMERIC,locale);
}



bool GaussianFileParser::ExtractJBlock(std::vector<QuantumCoupling>& c)
{
    //get the header
    std::string line;
    std::getline(*m_file,line);
    StringTokenizer tok(line," \t\r");
    std::vector<int> col(tok.size());
    for(size_t i=0;i<tok.size();++i)
    {
        col[i]=atoi(tok[i].c_str());
    }
    std::streampos pos;
    while(std::getline(*m_file,line) )
    {

        StringTokenizer tok(line," \t\r");
        if ( !kryomol::isnum(tok.front()) ) return false;
        if ( !kryomol::isinteger( tok.back()) )
        {
            int row=std::atoi(tok.front().c_str());
            for(size_t i=1;i<tok.size();++i)
                c.push_back(QuantumCoupling(row-1,col.at(i-1)-1,kryomol::atof(tok.at(i))));
        }
        else
        {
            m_file->seekg(pos);
            return true;
        }
        pos=m_file->tellg();
    }
    return false;
}



kryomol::D2Array<double> GaussianFileParser::ParseMagneticSusceptibility()
{
    const char* locale=std::setlocale(LC_NUMERIC,"C");
    m_file->clear();
    m_file->seekg(0,std::ios::beg);

    std::string line;
    kryomol::D2Array<double> t;
    bool bfound=false;
    while( std::getline(*m_file,line) )
    {
        if ( line.find("Magnetic susceptibility (cgs-ppm)") != std::string::npos )
            bfound=true;
        if ( bfound ) break;
    }
    if ( !bfound ) return t;
    std::getline ( *m_file,line );
    if ( line.find ( "Isotropic =" ) != std::string::npos )
    {
        t=D2Array<double>(3,3);

        for ( int i=0;i<3;++i )
        {
            std::getline ( *m_file,line );
            StringTokenizer token ( line," \t\r" );
            t ( 0,i ) =std::atof ( token.at ( 1 ).c_str() );
            t ( 1,i ) =std::atof ( token.at ( 3 ).c_str() );
            t ( 2,i ) = std::atof ( token.at ( 5 ).c_str() );
        }

    }
    std::setlocale(LC_NUMERIC,locale);
    return t;

}


bool GaussianFileParser::GetDipole ( const std::streampos& begin, const std::streampos& end )
{
    m_file->clear ( );               //
    m_file->seekg ( begin, std::ios::beg );   // moves pointer to

    std::string line;
    std::streampos pos=0;
    while ( std::getline ( *m_file,line ) )
    {
        StringTokenizer token ( line );
        if ( token.size() )
        {
            if ( token[0]=="Dipole" && token[1]=="moment" )
            {
                pos=m_file->tellg();
                goto lb1;
            }
        }
    }

lb1:

    if ( pos <= 0 ) return false;
    if ( end > 0 && pos > end ) return false;


    m_file->clear ( );               //

    m_file->seekg ( pos, std::ios::beg );   // moves pointer to

    std::getline ( *m_file,line );
    StringTokenizer token ( line );
    if ( token.size() < 6 ) return false; //simply return on error
    //change from physicist to chemist convention
    Coordinate c;
    c.x()=-std::atof ( token[1].c_str() );
    c.y()=-std::atof ( token[3].c_str() );
    c.z()=-std::atof ( token[5].c_str() );
    Molecules()->back().Frames().back().SetDipole(c);

    return true;

}
