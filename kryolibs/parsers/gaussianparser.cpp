/*****************************************************************************************
                            gaussianparser.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "gaussianparser.h"
#include "stringtools.h"
#include "molecule.h"
#include "exception.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <clocale>

using namespace kryomol;

GaussianParser::GaussianParser ( const char* file ) : Parser ( file )
{}

GaussianParser::GaussianParser ( std::istream* stream ) : Parser ( stream )
{}

GaussianParser::~GaussianParser()
{}

bool GaussianParser::ParseFile ( std::streampos pos )
{

    m_file->clear();
    m_file->seekg ( pos,std::ios::beg );
    m_file->clear();
    m_file->seekg ( pos,std::ios::beg );

    GetGeometry();
    m_file->clear();
    m_file->seekg ( pos,std::ios::beg );

    return true;

}

bool GaussianParser::GetGeometry()
{
    std::string line;
    std::vector<std::streampos> pos;
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

                pos.push_back ( m_file->tellg() );

            }

        }
        //Dont bother with jobs get the last
        //if( line.find('#',0) == 1 ) std::cout << "breaking" << std::endl;
        // if( line.find('#',0) == 1) break; //split different jobs

    }
    //Kein Standard Orientation, try input
    if ( pos.empty() )
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

                pos.push_back ( m_file->tellg() );

            }
            if ( !line.empty() )
                if ( line.at(0)  == '#'  ) break; //split different jobs
        }

    }
    if ( pos.empty() ) return false;

    std::vector<std::streampos>::iterator pt;
    //Get simlply the last conformer
    Molecules()->push_back ( Molecule() );
    Molecules()->back().Frames().push_back ( Frame ( &Molecules()->back() ) );
    for ( pt=pos.end()-1;pt!=pos.end();pt++ )
    {
        m_file->clear();
        m_file->seekg ( *pt,std::ios::beg );
        std::getline ( *m_file,line ); //skip one line
        while ( std::getline ( *m_file,line ) )
        {
            StringTokenizer token ( line," " );
            size_t size=token.size();

            if ( size >= 5 )
            {
                std::cout << "pusing back atom" << std::endl;
                Molecules()->back().Atoms().push_back ( Atom ( atoi ( token.at ( 1 ).c_str() ) ) );
                Coordinate c;
                c.x() =std::atof ( token.at ( size-3 ).c_str() );
                c.y() =std::atof ( token.at ( size-2 ).c_str() );
                c.z() =std::atof ( token.at ( size-1 ).c_str() );
                Molecules()->back().Frames().back().XYZ().push_back ( c );
            }

            else
                break;
        }


    }

    return true;
}

bool GaussianParser::HasKeyword ( std::string& line )
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


std::vector<JobHeader>& GaussianParser::Jobs()
{
    //lets get the lines beginning with a #
    std::string line;
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
        if ( line.find ( '#',0 ) == 1 ) //should be the second character
        {
            JobType tjob=singlepoint;
            bool bfound=false;

            std::string wholeline;
            while ( line.find ( "-----" ) == std::string::npos )
            {
                std::cout << "line now is " << line << std::endl;
                wholeline+=line;
                std::getline ( *m_file,line );
                line.erase ( 0,1 );
            }
            line=wholeline;
            //take into account that route can be splitted into several lines
            do
            {
                line=toupper ( line );
                if ( line.find ( "OPT",0 ) != std::string::npos )
                {
                    tjob=opt;
                    bfound=true;
                }
                else
                    if ( line.find ( "FREQ",0 ) != std::string::npos )
                    {
                        tjob=freq;
                        bfound=true;
                    }
                if ( line.find ( "TD",0 ) != std::string::npos )
                {
                    tjob=uv;
                    bfound=true;
                }
                if ( line.find ( "NMR",0 ) != std::string::npos )
                {
                    tjob= nmr;
                    bfound=true;
                }
                if ( line.find ( "ADMP",0 ) != std::string::npos )
                {
                    tjob= dyn;
                    bfound=true;
                }
                if ( bfound ) break;
                //pos=m_file->tellg();
                std::getline ( *m_file,line );
            }
            while ( line.find ( "-----" ) == std::string::npos );
            m_jobpos.push_back ( JobHeader( tjob,pos ) );
        }
        pos=m_file->tellg();
    }
    while ( std::getline ( *m_file,line ) );

    return m_jobpos;


}

void GaussianParser::ParseChemicalShifts()
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

void GaussianParser::ParseCouplingConstants(std::vector<QuantumCoupling>& c)
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


bool GaussianParser::ExtractJBlock(std::vector<QuantumCoupling>& c)
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

kryomol::D2Array<double> GaussianParser::ParseMagneticSusceptibility()
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
    std::cout << "in pms" << std::endl;
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
