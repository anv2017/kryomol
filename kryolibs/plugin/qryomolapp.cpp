/*****************************************************************************************
                            qryomolapp.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <QMessageBox>

#include "qryomolapp.h"
#include "pluginloader.h"
#include "applicationinfo.h"
#include <iostream>
#include <sstream>
#include <stdlib.h>



class kryomol::KryoMolApplicationPrivate
{
   public:
   KryoMolApplicationPrivate()
   {
     m_loader= new PluginLoader();
   }
   ~KryoMolApplicationPrivate()
    {
       delete m_loader;
    }
   kryomol::PluginLoader* m_loader;

};

QString m_file;

using namespace kryomol;

/** \brief KryoMol main application constructor

Bult a GUI enabled QryoMol main application
*/
KryoMolApplication::KryoMolApplication ( int & argc, char** argv ) : QApplication(argc,argv)
{
  _d = new KryoMolApplicationPrivate();
  
  //setOrganizationName(orgname);
  //setApplicationName(aname);
}

/** \brief KryoMol main application constructor

 In order to be used as a console application QryoMolApplication should be built with \a GUIenabled set to false
 Currently application can be only used in the graphical mode
*/
KryoMolApplication::KryoMolApplication( int & argc, char ** argv, bool GUIenabled ) : QApplication(argc,argv,GUIenabled)
{
 _d = new KryoMolApplicationPrivate();
  //setOrganizationName(orgname);
  //setApplicationName(aname);
}

KryoMolApplication::~KryoMolApplication()
{
  delete _d;
}

/** \return The local path of the currently loaded molecular structure file*/
QString KryoMolApplication::File()
{
  return m_file;
}

/** Set the path to the currently loaded molecular structure file*/
void KryoMolApplication::SetFile(const QString& file)
{
  m_file=file;
}

/** Set the path to the currently loaded molecular structure file*/
void KryoMolApplication::SetFile(const char* file)
{
  m_file=file;
}

/** 
 \return a stl vector with pointers to the found plugins
 */
std::vector<kryomol::Plugin*>& KryoMolApplication::Plugins()
{
  return _d->m_loader->Plugins();
}


/** \brief initialize the plugins the QryoMol plugins

This method will load the plugins at runtime using a PluginLoader object
*/
void KryoMolApplication::InitializePlugins()
{
  _d->m_loader->Load();
}

/** \brief Load the KryoMol plugins

This method will load the plugins at runtime using a PluginLoader object
*/
void KryoMolApplication::LoadPlugins()
{
  _d->m_loader->Load();
}

QString KryoMolApplication::Authors() const
{
  //return kryomol::authors;
}

QString KryoMolApplication::URL() const
{
  //return kryomol::url;
}
QString KryoMolApplication::Version() const
{
  //return kryomol::version;
}

bool KryoMolApplication::notify(QObject* receiver, QEvent* e)
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
