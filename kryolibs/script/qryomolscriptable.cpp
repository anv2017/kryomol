/*****************************************************************************************
                            qryomolscriptable.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "qryomolscriptable.h"
#include "stringtools.h"
#include "world.h"
#include "molecule.h"
#include "glvisor.h"
#include "exception.h"
#include "qryomolapp.h"

#include <QMessageBox> 

using namespace kryomol;

KryoMolScriptable::KryoMolScriptable()
{
}

void KryoMolScriptable::popmessage(QString text)
{
    QMessageBox k;
    k.setText(text);
    k.exec();
}



void KryoMolScriptable::superimpose(QString text)
{
    std::string str=text.toStdString();
    if ( str.empty() )
    {}
    else
    {
        StringTokenizer tok(str,",");
        std::vector<size_t> atomlist;
        StringTokenizer::iterator it;
        for(it=tok.begin();it!=tok.end();++it)
        {
            StringTokenizer tok1(*it,"-");
            if ( tok1.size() == 1 )
            {
                atomlist.push_back(atoi(tok1.at(0).c_str())-1);
            }
            else if ( tok1.size() == 2 )
            {
                int initial=atoi( tok1.at(0).c_str())-1;
                int final=atoi(tok1.at(1).c_str())-1;
                std::cout << "final is " << final << std::endl;
                for(int i=initial;i<=final;++i)
                {
                    atomlist.push_back(i);
                }
            }
            else
            {
                throw kryomol::Exception("Invalid syntaxis");
            }

        }
        std::cout << "atom list size is" << std::endl;
        if ( atomlist.size() < 3 )
            throw kryomol::Exception("Need at least three atoms to align frames");

        size_t refframe=m_world->CurrentMolecule()->CurrentFrameIndex();
        m_world->CurrentMolecule()->SuperImpose(refframe,atomlist);
        std::cout << "The following atoms have been superimposed:" << "(" << std::endl;
        for(std::vector<size_t>::const_iterator it=atomlist.begin();it!=atomlist.end();++it)
        {
            std::cout << "," << ((*it)+1);
        }
        std::cout << ")" << std::endl;

        m_app->InitializePlugins();
        m_world->Visor()->CenterMolecule(true,true);

    }
}

void KryoMolScriptable::eckarttransform(QString text)
{
    std::string str=text.toStdString();
    if ( str.empty() )
    {}
    else
    {
        StringTokenizer tok(str,",");
        std::vector<size_t> atomlist;
        StringTokenizer::iterator it;
        for(it=tok.begin();it!=tok.end();++it)
        {
            StringTokenizer tok1(*it,"-");
            if ( tok1.size() == 1 )
            {
                atomlist.push_back(atoi(tok1.at(0).c_str())-1);
            }
            else if ( tok1.size() == 2 )
            {
                int initial=atoi( tok1.at(0).c_str())-1;
                int final=atoi(tok1.at(1).c_str())-1;
                std::cout << "final is " << final << std::endl;
                for(int i=initial;i<=final;++i)
                {
                    atomlist.push_back(i);
                }
            }
            else
            {
                throw kryomol::Exception("Invalid syntaxis");
            }

        }
        std::cout << "atom list size is" << std::endl;
        if ( atomlist.size() < 3 )
            throw kryomol::Exception("Need at least three atoms to align frames");

        size_t refframe=m_world->CurrentMolecule()->CurrentFrameIndex();
        std::cout << "The following atoms have been superimposed:" << "(" << std::endl;
        for(std::vector<size_t>::const_iterator it=atomlist.begin();it!=atomlist.end();++it)
        {
            std::cout << "," << ((*it)+1);
        }
        m_world->CurrentMolecule()->EckartTransform(refframe,atomlist);
        std::cout << "hola mundo" << std::endl;
        std::cout << "The following atoms have been superimposed:" << "(" << std::endl;
        for(std::vector<size_t>::const_iterator it=atomlist.begin();it!=atomlist.end();++it)
        {
            std::cout << "," << ((*it)+1);
        }
        std::cout << ")" << std::endl;

        m_app->InitializePlugins();
        //m_world->Visor()->Initialize();
        m_world->Visor()->CenterMolecule(true);
        // m_world->Visor()->update();

    }
}

void KryoMolScriptable::SetWorld(World* w)
{
    m_world=w;
}

void KryoMolScriptable::SetApplication(KryoMolApplication* app)
{
    m_app=app;
}
