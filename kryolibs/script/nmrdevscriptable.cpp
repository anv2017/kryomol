/*****************************************************************************************
                            nmrdevscriptable.cpp  -  description
                             -------------------
    copyright            : (C) 2011 by Armando Navarro-Vazquez and Noa Campos-Lopez
    email                : armando.navarro@uvigo.es, noa.campos@uvigo.es
******************************************************************************************/

#include "nmrdevscriptable.h"
#include "stringtools.h"
#include "world.h"
#include "molecule.h"
#include "glvisor.h"
#include "exception.h"
#include "nmrdevapp.h"

#include <QMessageBox> 
using namespace nmrdev;

NMRDevScriptable::NMRDevScriptable()
{
}

void NMRDevScriptable::popmessage(QString text)
{
  QMessageBox k;
  k.setText(text);
  k.exec();
}



void NMRDevScriptable::superimpose(QString text)
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
            throw nmrdev::Exception("Invalid syntaxis");
          }

        }
        std::cout << "atom list size is" << std::endl;
        if ( atomlist.size() < 3 )
          throw nmrdev::Exception("Need at least three atoms to align frames");

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

void NMRDevScriptable::eckarttransform(QString text)
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
            throw nmrdev::Exception("Invalid syntaxis");
          }

        }
        std::cout << "atom list size is" << std::endl;
        if ( atomlist.size() < 3 )
          throw nmrdev::Exception("Need at least three atoms to align frames");

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

void NMRDevScriptable::SetWorld(World* w)
{
  m_world=w;
}

void NMRDevScriptable::SetApplication(NMRDevApplication* app)
{
  m_app=app;
}
