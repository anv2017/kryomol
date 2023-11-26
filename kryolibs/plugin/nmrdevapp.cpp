/*****************************************************************************************
                            nmrdevapp.cpp  -  description
                             -------------------
    copyright            : (C) 2011 by Armando Navarro-Vazquez and Noa Campos-Lopez
    email                : armando.navarro@uvigo.es, noa.campos@uvigo.es
******************************************************************************************/

#include <QMessageBox>

#include "nmrdevapp.h"
#include "pluginloader.h"
#include "applicationinfo.h"
#include <iostream>
#include <sstream>
#include <stdlib.h>



class nmrdev::NMRDevApplicationPrivate
{
   public:
   NMRDevApplicationPrivate()
   {
     m_loader= new PluginLoader();
   }
   ~NMRDevApplicationPrivate()
    {
       delete m_loader;
    };
   nmrdev::PluginLoader* m_loader;

};

QString m_file;
using namespace nmrdev;

/** \brief NMRDev main application constructor

Bult a GUI enabled NMRDev main application
*/
NMRDevApplication::NMRDevApplication ( int & argc, char** argv ) : QApplication(argc,argv)
{
  _d = new NMRDevApplicationPrivate();
  
  setOrganizationName(orgname);
  setApplicationName(aname);
}

/** \brief NMRDev main application constructor

 In order to be used as a console application NMRDevApplication should be built with \a GUIenabled set to false
 Currently application can be only used in the graphical mode
*/
NMRDevApplication::NMRDevApplication( int & argc, char ** argv, bool GUIenabled ) : QApplication(argc,argv,GUIenabled)
{
 _d = new NMRDevApplicationPrivate();
  setOrganizationName(orgname);
  setApplicationName(aname);
}

NMRDevApplication::~NMRDevApplication()
{
  delete _d;
}

/** \return The local path of the currently loaded molecular structure file*/
QString NMRDevApplication::File()
{
  return m_file;
}

/** Set the path to the currently loaded molecular structure file*/
void NMRDevApplication::SetFile(const QString& file)
{
  m_file=file;
}

/** Set the path to the currently loaded molecular structure file*/
void NMRDevApplication::SetFile(const char* file)
{
  m_file=file;
}

/** 
 \return a stl vector with pointers to the found plugins
 */
std::vector<nmrdev::Plugin*>& NMRDevApplication::Plugins()
{
  return _d->m_loader->Plugins();
}


/** \brief initialize the plugins the NMRDev plugins

This method will load the plugins at runtime using a PluginLoader object
*/
void NMRDevApplication::InitializePlugins()
{
  _d->m_loader->Load();
}

/** \brief Load the NMRDev plugins

This method will load the plugins at runtime using a PluginLoader object
*/
void NMRDevApplication::LoadPlugins()
{
  _d->m_loader->Load();
}

QString NMRDevApplication::Authors() const
{
  return nmrdev::authors;
}

QString NMRDevApplication::URL() const
{
  return nmrdev::url;
}
QString NMRDevApplication::Version() const
{
  return nmrdev::version;
}

bool NMRDevApplication::notify(QObject* receiver, QEvent* e)
{
  try
  {
    return QApplication::notify(receiver,e);
  }
  catch(std::exception& e)
  {
     QMessageBox k;
     k.setText(QString(e.what()));
     k.exec();
	 return true;
  }

}
