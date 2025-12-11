/*****************************************************************************************
                            orcaparser.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "orcaparser.h"
#include "stringtools.h"
#include "molecule.h"
#include "exception.h"
#include "stdlib.h"
#include "mathtools.h"
#include "orbitalarray.h"

#include <QProgressDialog>

#include <vector>
#include <iostream>
#include <cmath>
#include <clocale>

using namespace kryomol;

OrcaParser::OrcaParser ( const char* file ) : Parser ( file )
{

}

OrcaParser::OrcaParser ( std::istream* stream ) : Parser ( stream )
{
}

OrcaParser::~OrcaParser()
{}

bool OrcaParser::ParseFile ( std::streampos pos )
{
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

bool OrcaParser::GetGeometry()
{
    std::string line;
    bool bfound=false;
    while ( std::getline ( *m_file,line ) )
    {
        if ( line.find ( "CARTESIAN COORDINATES (ANGSTROEM)" ) != std::string::npos )
        {
            bfound=true;

            std::getline ( *m_file,line );

            m_pos.push_back ( m_file->tellg() );
        }
        if ( line.find("JOB NUMBER") != std::string::npos ) break;
    }
    if ( m_pos.empty() ) return false;

    std::vector<std::streampos>::iterator pt;

    //Get simply the last molecule
    Molecules()->push_back ( Molecule() );
    pt=m_pos.begin();
    m_file->clear();
    m_file->seekg ( *pt,std::ios::beg );
    while ( std::getline ( *m_file,line ) )
    {
        StringTokenizer token ( line," " );
        size_t size=token.size();

        if ( size >= 4 )
        {
            Molecules()->back().Atoms().push_back ( Atom ( token.at ( 0 ).c_str() ) );
        }
        else break;
    }
    for ( pt=m_pos.begin();pt!=m_pos.end();pt++ )
    {
        Molecules()->back().Frames().push_back ( Frame ( &Molecules()->back() ) );

        m_file->clear();
        m_file->seekg ( *pt,std::ios::beg );
        while ( std::getline ( *m_file,line ) )
        {
            StringTokenizer token ( line," " );
            size_t size=token.size();

            if ( size >= 4 )
            {
                Coordinate c;
                c.x() =std::atof ( token.at ( size-3 ).c_str() );
                c.y() =std::atof ( token.at ( size-2 ).c_str() );
                c.z() =std::atof ( token.at ( size-1 ).c_str() );
                Molecules()->back().Frames().back().XYZ().push_back ( c );
            }
            else  break;
        }
        GetEnergyForFrame();
        GetGradientForFrame();

        //Get the SCF energy

    }

    Molecules()->back().SetBonds();

    //just keepe the last position in case we needed later
    //since Orca do not write agian geoemtyr in the job section
    if ( m_pos.size() > 1 )
    {
        m_pos.erase(m_pos.begin(), m_pos.end() - 1);
    }


    return true;
}

void OrcaParser::GetEnergyForFrame()
{
    std::string line;
    while(std::getline(*m_file,line) )
    {
        if ( line.find( "FINAL SINGLE POINT ENERGY") != std::string::npos )
        {
            StringTokenizer tok(line," \t:");
            if ( tok.size() < 3 ) throw kryomol::Exception("error parsing SCF energy");
            Molecules()->back().Frames().back().SetEnergy(std::stod(tok.back()),"SCF");
            return;
        }
    }
}

void OrcaParser::GetGradientForFrame()
{
    std::string line;
    while(std::getline(*m_file,line) )
    {
        if ( line.find( "|Geometry convergence|") != std::string::npos )
        {
            while(std::getline(*m_file,line))
            {
                if ( line.find("RMS gradient") != std::string::npos )
                {
                    StringTokenizer tok(line," \t:");
                    Molecules()->back().Frames().back().SetRMSForce(std::stod(tok.at(2)));

                }

                if ( line.find("MAX gradient") != std::string::npos )
                {
                    StringTokenizer tok(line," \t:");
                    Molecules()->back().Frames().back().SetMaximumForce(std::stod(tok.at(2)));

                }

                if ( line.find("RMS step") != std::string::npos )
                {
                    StringTokenizer tok(line," \t:");
                    Molecules()->back().Frames().back().SetRMSDisplacement(std::stod(tok.at(2)));

                }

                if ( line.find("MAX step") != std::string::npos )
                {
                    StringTokenizer tok(line," \t:");
                    Molecules()->back().Frames().back().SetMaximumDisplacement(std::stod(tok.at(2)));
                    return;

                }

            }


        }
    }
}


bool OrcaParser::ParseOrbitals(std::streampos pos)
{
    bool b = false;

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

            b = GetOrbitalData();
        }
    }

    return b;
}

bool OrcaParser::GetBasisCenters()
{
    std::string line;
    std::vector<size_t> v_atoms;
    std::vector< std::vector<Orbital> > v_orbitals;
    while ( std::getline ( *m_file,line ) )
    {
        m_norbitals = 0;
        if ( line.find ( "BASIS SET IN INPUT FORMAT" ) != std::string::npos )
        {
            std::getline( *m_file,line );
            std::getline( *m_file,line );
            std::getline( *m_file,line );
            while (line.find ( "# Basis set for element :" ) != std::string::npos)
            {
                int norbitals = 0;
                StringTokenizer tokenAtom (line, " \t\r");
                std::string atom = tokenAtom.at(tokenAtom.size()-1);
                std::vector<Orbital> orbitals;

                std::getline( *m_file,line );
                std::getline( *m_file,line );

                StringTokenizer token ( line," \t\r" );
                size_t sizeOrbital = token.size();
                while (sizeOrbital == 2)
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
                                xs.push_back(atof(tokenCoefficient.at(2)));
                                xp.push_back(atof(tokenCoefficient.at(3)));
                                alpha.push_back(atof(tokenCoefficient.at(1)));
                            }
                        }
                        else
                        {
                            if (tokenCoefficient.size()>=2)
                            {
                                xs.push_back(atof(tokenCoefficient.at(2)));
                                alpha.push_back(atof(tokenCoefficient.at(1)));
                            }
                        }
                    }

                    if (tokenOrbital.at(0)=="S")
                    {
                        orbitals.push_back( Orbital (Orbital::S,alpha,xs,xp));
                        norbitals += 1;
                    }
                    if (tokenOrbital.at(0)=="SP")
                    {
                        orbitals.push_back( Orbital (Orbital::SP,alpha,xs,xp));
                        norbitals += 4;
                    }
                    if (tokenOrbital.at(0)=="P")
                    {
                        orbitals.push_back(Orbital (Orbital::P,alpha,xs,xp));
                        norbitals += 3;
                    }
                    if (tokenOrbital.at(0)=="D")
                    {
                        orbitals.push_back(Orbital (Orbital::D,alpha,xs,xp));
                        norbitals += 5;
                    }
                    if (tokenOrbital.at(0)=="F")
                    {
                        orbitals.push_back(Orbital (Orbital::F,alpha,xs,xp));
                        norbitals += 7;
                    }

                    std::getline( *m_file,line );
                    StringTokenizer token ( line," \t\r" );
                    sizeOrbital = token.size();
                }

                std::getline( *m_file,line );
                std::getline( *m_file,line );

                std::vector<Atom>::iterator it;
                size_t at = 0;
                for (it = Molecules()->back().Atoms().begin(); it != Molecules()->back().Atoms().end(); ++it, ++at)
                {
                    if (it->Symbol() == atom)
                    {
                        m_norbitals += norbitals;
                        v_atoms.push_back(at);
                        v_orbitals.push_back(orbitals);
                    }
                }

            }
            for (size_t a=0; a<v_atoms.size(); ++a)
            {
                size_t it=0;
                while ((v_atoms.at(it)!=a) && (it<v_atoms.size()))
                {
                    ++it;
                }
                m_atoms.push_back(a);
                m_orbitals.push_back(v_orbitals.at(it));
            }

            return true;
        }
    }
    return false;
}

bool OrcaParser::GetOrbitalData()
{
    bool b =false;
    int frame =0;
    std::vector<Frame>::iterator ft = Molecules()->back().Frames().begin();
    size_t pt=0;
    std::string line;
    while ( std::getline(*m_file,line) )
    {
        if ( line.find ( "MOLECULAR ORBITALS" ) != std::string::npos )
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
            std::getline( *m_file,line );
            int K=0;
            while ( !(line.empty()) )
            {
                StringTokenizer tokenSize ( line," " );
                int M=tokenSize.size();

                std::getline( *m_file,line );
                StringTokenizer tokenEigenv ( line," " );
                for (int e=0; e<M; e++)
                {
                    eigenvalues.push_back(atof(tokenEigenv.at(tokenEigenv.size()-M+e)));
                }
                std::getline( *m_file,line );
                StringTokenizer tokenOcc ( line," " );
                for (int e=0; e<M; e++)
                {
                    occupations.push_back(atof(tokenOcc.at(tokenOcc.size()-M+e)));
                }
                std::getline( *m_file,line );
                for (int n=0; n<m_norbitals; n++)
                {
                    std::getline( *m_file,line );
                    StringTokenizer tokenCoeff ( line," " );
                    for (int m=0; m<M; m++)
                    {
                        matrix(n,K+m) = atof(tokenCoeff.at(tokenCoeff.size()-M+m));
                        //qDebug() << "(" << n << "," << K+m << ") = " << atof(tokenCoeff.at(tokenCoeff.size()-M+m)) << endl;
                    }
                }
                K += M;
                std::getline( *m_file,line );
            }

            //Reorder the molecular orbitals coefficients
            size_t orbital=0;
            for (std::vector< std::vector<Orbital> >::iterator it=m_orbitals.begin(); it!=m_orbitals.end(); ++it)
            {
                for (std::vector<Orbital>::iterator ot=(*it).begin(); ot!=(*it).end(); ++ot)
                {
                    switch ((*ot).Type())
                    {
                        qDebug() << "Type: " << (*ot).Type() << endl;
                    case Orbital::S:
                        orbital = orbital+1;
                        break;
                    case Orbital::SP:
                        matrix.SwapRows(orbital+1,orbital+2);
                        matrix.SwapRows(orbital+2,orbital+3);
                        orbital = orbital+4;
                        break;
                    case Orbital::P:
                        matrix.SwapRows(orbital,orbital+1);
                        matrix.SwapRows(orbital+1,orbital+2);
                        orbital = orbital+3;
                        break;
                    case Orbital::D:
                        orbital = orbital+5;
                        break;
                    case Orbital::F:
                        orbital = orbital+7;
                        break;
                    }
                }
            }

            std::vector<BasisCenter> basis;
            for (size_t i=0; i<m_atoms.size(); ++i)
                basis.push_back(BasisCenter(ft->XYZ().at(m_atoms.at(i)),m_orbitals.at(i)));
            ft->OrbitalsData().SetHomo(m_homo);
            ft->OrbitalsData().SetLumo(m_lumo);
            ft->OrbitalsData().SetTypeD(5);
            ft->OrbitalsData().SetTypeF(7);
            ft->OrbitalsData().SetBasisCenters(basis);
            ft->OrbitalsData().SetCoefficients(matrix);
            ft->OrbitalsData().SetEigenvalues(eigenvalues);
            b = true;

            ++pt;
            ++ft;
            frame+=1;
        }
    }
    return b;
}

bool OrcaParser::ExistOrbitals()
{
    std::string line;
    while ( std::getline ( *m_file,line ))
    {
        if ( line.find("MOLECULAR ORBITALS") != std::string::npos )
            return true;
    }
    return false;
}

bool OrcaParser::GetHomoLumo()
{
    std::string line;
    while ( std::getline ( *m_file,line ) )
    {
        if ( line.find ( "Number of Electrons" ) != std::string::npos )
        {
            StringTokenizer token(line," \t\r");
            m_homo = floor(atof(token.at(5))/2);
            m_lumo = m_homo+1;

            return true;
        }
    }
    return false;
}

std::vector<JobHeader>& OrcaParser::Jobs()
{
    std::string line;
    m_file->clear();
    m_file->seekg(0);
    std::streampos pos=m_file->tellg();


    while(std::getline(*m_file,line) )
    {
        if ( line.find("* Geometry Optimization Run *") != std::string::npos
            || line.find("*    Relaxed Surface Scan    *") != std::string::npos )
        {
            pos=m_file->tellg();
            m_jobpos.push_back(JobHeader(opt,pos));
        }

        //if ( line.find("*     ORCA property calculations      *") != std::string::npos )
        if ( line.find("* Single Point Calculation *") != std::string::npos ||
            line.find("*     ORCA property calculations      *") != std::string::npos ||
            line.find("Energy+Gradient Calculation")!= std::string::npos )
        {
            pos=m_file->tellg();
            m_jobpos.push_back(JobHeader(singlepoint,pos));
            break;
        }
    }

    //just one job
    if ( m_jobpos.empty() )
    {
        m_jobpos.push_back(JobHeader(singlepoint,0));
    }

    for(std::vector<JobHeader>::iterator it=m_jobpos.begin();it!=m_jobpos.end();++it)
    {
        if ( it->type == singlepoint )
        {
            m_file->clear();
            m_file->seekg(it->pos);
            while(std::getline(*m_file,line) )
            {
                bool bfound=false;


                if ( line.find("TD-DFT/TDA EXCITED STATES") != std::string::npos
                    || line.find("TD-DFT EXCITED STATES") != std::string::npos  )
                {
                    it->type=uv;
                    bfound=true;
                }

                if ( line.find("ORCA NUMERICAL FREQUENCIES") != std::string::npos || line.find("ORCA SCF HESSIAN") != std::string::npos )
                {
                    it->type=freq;
                    bfound=true;
                }
                if ( bfound ) break;
            }
        }

    }

    return m_jobpos;
}

bool OrcaParser::ParseFrequencies ( std::streampos pos )
{
    m_file->clear();
    m_file->seekg ( pos,std::ios::beg );
    if ( !GetFrequencies() ) return false;
    m_file->clear();
    m_file->seekg ( pos,std::ios::beg );

    Molecule& molecule=Molecules()->back();
    molecule.Frames().back().AllocateVectors(true);

    GetNormalModes();

    /*m_file->clear();
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
    qDebug() << "LIST FREQUENCIES in ParseFrequencies: " << v.size() << endl;*/

    return true;
}

bool OrcaParser::GetFrequencies()
{
    Molecule& molecule=Molecules()->back();

    std::vector<Frequency>& frequencies=molecule.Frames().back().GetFrequencies();

    std::string line;

    while ( std::getline ( *m_file,line ) )
    {
        if ( line.find("Scaling factor for frequencies") != std::string::npos ) break;
    }

    std::getline ( *m_file,line );

    while( std::getline(*m_file,line) )
    {
        StringTokenizer tok(line," \t\r");
        if ( tok.empty() ) break;
        frequencies.push_back(Frequency(std::stof(tok.at(1)),0,0));
    }

    bool irfound=false;
    while( std::getline(*m_file,line) )
    {
        if (line.find("IR SPECTRUM") != std::string::npos ) irfound=true;
        if (irfound) break;
    }
    if ( irfound ==false ) return false;
    //skip the "---" closign line
    std::getline(*m_file,line);
    size_t Intcolumn=0;
    do
    {
        std::getline(*m_file,line);
        if ( line.find( "Mode") != std::string::npos )
        {
            StringTokenizer tok(line," \t\r");
            for(;Intcolumn<tok.size();++Intcolumn)
            {
                if ( tok[Intcolumn] == "Int" ) break;
            }
        }
    } while( line.find("---") == std::string::npos );

    while( std::getline(*m_file,line) )
    {
        StringTokenizer tok(line," \t\r");
        if ( tok.empty() ) break;
        //Intensity in km/mol
        frequencies.at(std::stoi(tok.at(0))).y=std::stof(tok.at(Intcolumn));
        //std::cout << "frequencies.y" << frequencies.at(std::stoi(tok.at(0))).y << std::endl;
    }

    //########
    bool ramanfound=false;
    while( std::getline(*m_file,line) )
    {
        if (line.find("RAMAN SPECTRUM") != std::string::npos ) ramanfound=true;
        if (ramanfound) break;
    }
    if ( ramanfound == false ) return false;
    //skip the "---" closign line
    std::getline(*m_file,line);
    Intcolumn=0;
    do
    {
        std::getline(*m_file,line);
        if ( line.find( "Mode") != std::string::npos )
        {
            StringTokenizer tok(line," \t\r");
            for(;Intcolumn<tok.size();++Intcolumn)
            {
                if ( tok[Intcolumn] == "Activity" ) break;
            }
        }
    } while( line.find("---") == std::string::npos );

    while( std::getline(*m_file,line) )
    {
        StringTokenizer tok(line," \t\r");
        if ( tok.empty() ) break;
        //Intensity in km/mol
        frequencies.at(std::stoi(tok.at(0))).w=std::stof(tok.at(Intcolumn));
        //std::cout << "frequencies.y" << frequencies.at(std::stoi(tok.at(0))).y << std::endl;
    }

    return true;

}

bool OrcaParser::GetNormalModes()
{
    bool modesfound=false;
    std::string line;
    while(std::getline(*m_file,line))
    {
        if (line.find("NORMAL MODES") != std::string::npos ) modesfound=true;
        if ( modesfound == true ) break;
    }
    if ( modesfound == false ) return false;

    while(std::getline(*m_file,line))
    {
        StringTokenizer tok(line," \t\r");
        if ( !tok.empty() )
        {
            if ( kryomol::isinteger(tok.front() )  )
            {
                if ( GetNormalModeBlock(tok) ) break;
            }
        }
    }

    return true;

}

bool OrcaParser::GetNormalModeBlock(StringTokenizer& headertok)
{
    Molecule& molecule=Molecules()->back();
    std::vector< std::vector<Coordinate> >& modes = molecule.GetMode();
    std::string line;
    size_t firstmode=std::stoi(headertok.front());
    size_t lastmode=std::stoi(headertok.back());
    std::streampos pos=m_file->tellg();
    while(std::getline(*m_file,line))
    {
        StringTokenizer tok(line," \t\r");
        if ( tok.empty() )
        {
            m_file->clear();
            m_file->seekg(pos);
            return true;
        }
        if ( kryomol::isinteger(tok.back()))
        {
            m_file->clear();
            m_file->seekg(pos);
            return false;
        }

        pos=m_file->tellg();
        int cnumber=std::stoi(tok.front());
        int cidx=cnumber/3;
        int cmod=cnumber % 3;
        size_t tidx=0;
        for(size_t i=firstmode;i<=lastmode;++i)
        {

            D1Array<float>& dc=modes.at(i).at(cidx);
            dc(cmod)=std::stof(tok.at(++tidx));
        }

    }
    return false;
}


void OrcaParser::ParseUVLengthBlock(std::vector<Spectralline>& lines)
{
    std::string line;
    while(std::getline(*m_file,line))
    {
        StringTokenizer token(line," \t");
        if ( token.size() >= 8 )
        {

            Spectralline sline;
            sline.x=std::stof( *(token.rbegin()+5) );
            sline.y0=std::stof( *(token.rbegin()+4) );
            lines.push_back(sline);
        }
        else break;
    }
}

void OrcaParser::ParseUVVelocityBlock(std::vector<Spectralline>& lines)
{
    std::string line;
    std::vector<Spectralline>::iterator spt=lines.begin();
    while(std::getline(*m_file,line))
    {
        StringTokenizer token(line," \t");
        if ( token.size() >= 7 )
        {

            spt->SetRotatoryStrengthVelocity(std::stof(*(token.rbegin()+3)));
            ++spt;
        }
        else break;
    }
}

void OrcaParser::ParseCDLengthBlock(std::vector<Spectralline>& lines)
{
    std::string line;
    std::vector<Spectralline>::iterator spt=lines.begin();
    while(std::getline(*m_file,line))
    {
        StringTokenizer token(line," \t");
        if ( token.size() >= 7 )
        {
            spt->SetRotatoryStrengthLength(std::stof(*(token.rbegin()+3)));
            ++spt;
        }
        else break;
    }

}

void OrcaParser::ParseCDVelocityBlock(std::vector<Spectralline>& lines)
{
    std::string line;
    std::vector<Spectralline>::iterator spt=lines.begin();
    while(std::getline(*m_file,line))
    {
        StringTokenizer token(line," \t");
        if ( token.size() >= 7 )
        {
            spt->SetRotatoryStrengthVelocity(std::stof(*(token.rbegin()+3)));
            ++spt;
        }
        else break;
    }

}

void OrcaParser::ParseSolventShiftBlock(std::vector<Spectralline>& lines)
{
    std::string line;
    std::vector<Spectralline>::iterator spt=lines.begin();
    while(std::getline(*m_file,line))
    {
        if ( line.find("State") != std::string::npos )
            break;
    }
    std::getline(*m_file,line);
    for(std::vector<Spectralline>::iterator spt=lines.begin();spt!=lines.end();++spt)
    {
        std::getline(*m_file,line);
        StringTokenizer token(line," \t");
        spt->SetSolventShift(std::atof(token.at(4).c_str()));
    }
}


bool OrcaParser::ParseUV ( std::streampos pos )
{
    m_file->clear();
    m_file->seekg ( pos,std::ios::beg );
    std::string line;
    Molecule& molecule=Molecules()->back();
    std::vector<Spectralline>& lines=molecule.Frames().back().GetSpectralLines();

    //There is a bug in Orca6 ? so only CD SPECTRUM appears for the lenght formalism
    bool cdlengthdone=false;
    while ( std::getline ( *m_file,line ) )
    {
        if ( line.find("ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS") != std::string::npos )
        {
            //skip four lines
            for (int i=0;i<4;++i)
            {
                std::getline(*m_file,line);
            }
            ParseUVLengthBlock(lines);
        }

        /* if ( line.find("ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS") != std::string::npos )
        {
            //skip four lines
            for (int i=0;i<4;++i)
            {
                std::getline(*m_file,line);
            }
            ParseUVelocityBlock(lines);
        }*/

        if ( line.find("CD SPECTRUM") != std::string::npos && cdlengthdone == false)
        {
            //skip four lines
            for (int i=0;i<4;++i)
            {
                std::getline(*m_file,line);
            }
            ParseCDLengthBlock(lines);
            cdlengthdone=true;
        }


        if ( line.find("CD SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS") != std::string::npos )
        {
            //skip four lines
            for (int i=0;i<4;++i)
            {
                std::getline(*m_file,line);
            }
            ParseCDVelocityBlock(lines);
        }


    }

    m_file->clear();
    m_file->seekg ( pos,std::ios::beg );
    //I think this was probably a Orca4 thing for TD-DFT and it is no longer present in Orca5 or Orca6 save for %mdci computations
    while(std::getline(*m_file,line))
    {
        if ( line.find("CALCULATED SOLVENT SHIFTS") != std::string::npos )
        {
            //skip four lines
            for (int i=0;i<4;++i)
            {
                std::getline(*m_file,line);
            }
            ParseSolventShiftBlock(lines);
        }
    }
    //Get transition vectors
    //Get transition changes
    return true;
}

