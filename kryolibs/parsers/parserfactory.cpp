/*****************************************************************************************
                            parserfactory.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "parserfactory.h"
#include "parsers.h"
#include "stringtools.h"

using namespace kryomol;

ParserFactory::ParserFactory(const char* file)
{
    m_stream = new std::ifstream ( file,std::ios::binary ); m_bstreamcreated=true;
}

ParserFactory::ParserFactory(std::istream* stream) : m_stream(stream), m_bstreamcreated(false)
{

}

#ifdef __MINGW32__
ParserFactory::ParserFactory(std::filesystem::path p)
{
    m_stream = new std::ifstream ( p );
}
#endif


ParserFactory::~ParserFactory()
{
    if ( m_bstreamcreated) delete m_stream;
}

Parser* ParserFactory::BuildParser()
{
    Parser* p;
    switch ( this->GetFileType() )
    {
    case Aces:
        p= new AcesParser ( m_stream );
        break;
    case Gaussian:
        p= new GaussianParser ( m_stream );
        break;
    case GaussianArchive:
        p= new ArchiveParser ( m_stream );
        break;
    case GaussianInput:
        p = new GaussianInputParser ( m_stream );
        break;
    case GaussianFile:
        p = new GaussianFileParser (m_stream);
        break;
    case GaussianCube:
        p= new GaussianCubeParser ( m_stream );
        break;
    case Gamess:
        p= new GamessParser(m_stream );
        break;
    case Maestro:
        p = new MaestroParser ( m_stream );
        break;
    case MdlV2000:
        p= new MdlV2000Parser ( m_stream );
        break;
    case NwChem:
        p= new NwChemParser(m_stream );
        break;
    case PDB:
        p = new PDBParser ( m_stream );
        break;
    case MacroModel:
        p = new MacroModelParser ( m_stream );
        break;
    case XYZ:
        p = new XYZParser ( m_stream );
        break;
    case HyperChem:
        p = new HyperChemParser ( m_stream );
        break;
    case PCModel:
        p = new PCModelParser ( m_stream );
        break;
    case Orca:
        p = new OrcaParser ( m_stream );
        break;
    case None:
    default:
        p = NULL;
        break;
    }
    return p;
}

ParserFactory::filetype ParserFactory::GetFileType()
{

    //first try gaussian
    std::string line;

    m_stream->clear();
    m_stream->seekg(0,std::ios::beg);
    while(std::getline(*m_stream,line))
    {
        if(line.find("Standard orientation") != std::string::npos ||  line.find("Input orientation") != std::string::npos || line.find("Z-Matrix orientation") != std::string::npos )
        {
            return GaussianFile;
        }
    }

    m_stream->clear();
    m_stream->seekg(0,std::ios::beg);
    while(std::getline(*m_stream,line))
    {
        if(line.find("* O   R   C   A *") != std::string::npos )
        {
            return Orca;
        }
    }

    m_stream->clear();
    m_stream->seekg(0,std::ios::beg);
    while(std::getline(*m_stream,line))
    {
        if(line.find("GAMESS VERSION") != std::string::npos)
            return Gamess;
    }

    m_stream->clear();
    m_stream->seekg(0,std::ios::beg);
    while(std::getline(*m_stream,line))
    {
        if(line.find("s_m_m2io_version") != std::string::npos)
        {
            std::cout << "MacroModel Maestro format" << std::endl;
            return Maestro;
        }

    }

    m_stream->clear();
    m_stream->seekg ( 0,std::ios::beg );
    while ( std::getline ( *m_stream,line ) )
    {
        if ( line.find ( "ACES2: Advanced Concepts in Electronic Structure II" ) != std::string::npos || line.find("CFOUR Coupled-Cluster techniques for Computational Chemistry") != std::string::npos )
        {
            std::cout << "ACESII/CFOUR file" << std::endl;
            return Aces;
        }
    }


    m_stream->clear();
    m_stream->seekg(0,std::ios::beg);
    while ( std::getline ( *m_stream,line ) )
    {
        if ( line.find ( "Northwest Computational Chemistry Package" ) != std::string::npos )
        {
            std::cout << "NwChem file" << std::endl;
            return NwChem;
        }
    }


    m_stream->clear();
    m_stream->seekg(0,std::ios::beg);
    while ( std::getline ( *m_stream,line ) )
    {
        if ( line.find ( "[HIN System Description]" ) != std::string::npos )
        {
            std::cout << "HyperChem file" << std::endl;
            return HyperChem;
        }
    }


    m_stream->clear();
    m_stream->seekg(0,std::ios::beg);
    while ( std::getline ( *m_stream,line ) )
    {
        if ( line.find ( "{PCM " ) != std::string::npos )
        {
            std::cout << "PCModel file" << std::endl;
            return PCModel;
        }
    }


    m_stream->clear();
    m_stream->seekg(0,std::ios::beg);
    while ( std::getline ( *m_stream,line ) )
    {
        if ( line.find ( "MO coefficients" ) != std::string::npos )
        {
            std::cout << "GaussianCube file" << std::endl;
            return GaussianCube;
        }
    }


    m_stream->clear();
    m_stream->seekg(0,std::ios::beg);
    std::getline(*m_stream,line);
    StringTokenizer tok(line," \t\r,");
    if( tok.size()== 2 )
    {
        if ( isinteger(tok.at(0)) && isinteger(tok.at(1)) )
            return GaussianInput;
    }

    m_stream->clear();
    m_stream->seekg(0,std::ios::beg);

    //well, I dont know exactly what to search for, maybe 5 integers?
    while(std::getline(*m_stream,line) )
    {
        StringTokenizer token(line," \t\r");
        if ( token.size() > 4  )
        {
            std::cout << token.back() << std::endl;

            if  ( kryomol::toupper(token.back())  == "V2000" )
            {
                std::cout << "MDL V2000 file" << std::endl;
                return MdlV2000;
            }

        }
    }

    m_stream->clear();
    m_stream->seekg(0,std::ios::beg);
    std::getline(*m_stream,line);
    if( ( line.find( "1\\1\\" ) != std::string::npos  ) ||  ( line.find( "1|1|") != std::string::npos ) ) return GaussianArchive;

    m_stream->clear();
    m_stream->seekg(0,std::ios::beg);
    while(std::getline(*m_stream,line) )
    {
        if ( (line.find( "HETATM" ) != std::string::npos ) || line.find( "ATOM" ) != std::string::npos )
            return PDB;
    }

    m_stream->clear();
    m_stream->seekg(0,std::ios::beg);
    std::getline(*m_stream,line);
    std::getline(*m_stream,line);
    {
        StringTokenizer token(line," \t\r");
        if (token.size() > 13 )
        {
            for(int i=0;i<13;i++)
            {
                if ( !isinteger(token.at(i))) goto lb;
            }
            std::cout << "MacroModel Old Format" << std::endl;
            return MacroModel;
        }
    }
lb:
    m_stream->clear();
    m_stream->seekg(0,std::ios::beg);
    std::getline(*m_stream,line);
    StringTokenizer token(line," \t\r");
    if( token.size() == 1 )
    {
        if( kryomol::isnum(token.at(0)) ) return XYZ;
    }

    return None;

}

bool ParserFactory::isGaussianFile()
{
    bool ftype;

    if ( (this->GetFileType() == GaussianFile) || (this->GetFileType() == Aces) || ( this->GetFileType() == Orca ))
        ftype = true;
    else
        ftype = false;

    return ftype;
}


bool ParserFactory::existDensity()
{
    bool ftype = false;

    if (this->GetFileType() == GaussianCube)
        ftype = true;
    if (this->GetFileType() == GaussianFile)
    {
        std::string line;
        while ((std::getline(*m_stream,line)) && (!(ftype)))
        {
            if (line.find("AO basis set in the form of general basis input") != std::string::npos)
            {
                ftype = true;
            }
        }
    }
    if (this->GetFileType() == Orca)
    {
        std::string line;
        while ((std::getline(*m_stream,line)) && (!(ftype)))
        {
            if (line.find("BASIS SET IN INPUT FORMAT") != std::string::npos)
            {
                ftype = true;
            }
        }
    }
    return ftype;
}


bool ParserFactory::existOrbitals()
{
    bool ftype = false;

    if (this->GetFileType() == GaussianFile)
    {
        std::string line;
        while ((std::getline(*m_stream,line)) && (!(ftype)))
        {
            if (line.find("Molecular Orbital Coefficients" ) != std::string::npos)
            {
                ftype = true;
            }
        }
    }
    if (this->GetFileType() == Orca)
    {
        std::string line;
        while ((std::getline(*m_stream,line)) && (!(ftype)))
        {
            if (line.find("MOLECULAR ORBITALS") != std::string::npos)
            {
                ftype = true;
            }
        }
    }

    return ftype;
}

bool ParserFactory::existAlphaBeta()
{
    bool ftype = false;

    if (this->GetFileType() == GaussianFile)
    {
        std::string line;
        while ((std::getline(*m_stream,line)) && (!(ftype)))
        {
            if (line.find("Alpha Molecular Orbital Coefficients:" ) != std::string::npos)
            {
                ftype = true;
            }
        }
    }

    return ftype;
}
