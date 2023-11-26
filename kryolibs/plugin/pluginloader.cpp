/*****************************************************************************************
                            pluginloader.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <iostream>

#include <QDir>
#include <QLibrary>
#include <QCoreApplication>

#include "pluginloader.h"
#include "stringtools.h"
#include "plugin.h"



class kryomol::PluginLoaderPrivate
{
  public:
    friend class PluginLoader;
    PluginLoaderPrivate() {}
    ~PluginLoaderPrivate() {}
  private:
    QStringList m_dirs;
    std::vector<kryomol::Plugin*> m_plugins;
};

using namespace kryomol;

/** \brief PluginLoader constructor
    The list of directories containing qryomol plugins is filled inside
    this constuctor
*/
PluginLoader::PluginLoader()
{
  m_private = new PluginLoaderPrivate();
  InitializePluginDirs();
}


PluginLoader::~PluginLoader()
{
  delete m_private;
}

/** \return a stl vector with pointers to the succefully loaded plugins*/
std::vector<kryomol::Plugin*>& PluginLoader::Plugins()
{
  return m_private->m_plugins;
}

/**
 \return a list of the directories where application search for plugins
*/
const QStringList& PluginLoader::Dirs() const
{
  return m_private->m_dirs;
}

/** \brief Initialize the plugin directories

This function is called inside the constructor
*/
void PluginLoader::InitializePluginDirs()
{
  QString spldir;

#ifdef Q_WS_WIN

  QString baseplugindir;
  spldir = QCoreApplication::instance()->applicationDirPath();
  spldir.truncate ( spldir.lastIndexOf ( "bin" ) );
  spldir +="plugins";
  baseplugindir=spldir;
  m_private->m_dirs << baseplugindir;
#endif

#ifdef Q_OS_UNIX
  char* nppath=getenv ( "KRYOMOL_DIR" );

  QString baseplugindir ( nppath );
#ifdef Q_WS_X11
  baseplugindir+="/plugins";
#endif
#ifdef Q_OS_MAC
  QDir dir ( QCoreApplication::instance()->applicationDirPath() );
  dir.cdUp();
  dir.cd ( "PlugIns" );
  baseplugindir = dir.absolutePath();
#endif

  m_private->m_dirs << baseplugindir;
  nppath=getenv ( "KRYOMOL_PLUGINS_DIRS" );
  if ( nppath )
  {
    std::string paths ( nppath );
    std::cout << paths << std::endl;
    StringTokenizer token ( paths,":" );
    StringTokenizer::iterator it;

    for ( it=token.begin();it!=token.end();++it )
    {
      m_private->m_dirs << QString ( it->c_str() );
    }
  }

  QStringList::Iterator it;
  for ( it=m_private->m_dirs.begin();it!=m_private->m_dirs.end();++it )
  {
#ifdef Q_OS_MAC
    QCoreApplication::instance()->addLibraryPath ( *it );
#endif

  }

#endif
}

/** \brief load the qryomol plugins.

Search on predefined directories and load the plugins.
In order the plugin to be loaded, the GetPlugin() function should be succesfully resolved for the 
corresponding shared library*/
void PluginLoader::Load()
{
  QStringList::Iterator it;
  QString spldir;
  QDir pluginsdir;
  bool pluginfound=false;
#ifdef Q_OS_MAC
// m_private->m_dirs.clear();
// m_private->m_dirs << "/Users/armando/qryomol4/bin/qryomol.app/Contents/PlugIns";
#endif
  for ( it=m_private->m_dirs.begin();it!=m_private->m_dirs.end();++it )
  {
    std::cout << "Searching plugins on " << ( *it ).toStdString() << std::endl;


#ifdef Q_WS_WIN
    pluginsdir.setPath ( *it );
    QStringList namefilters;
    namefilters << "*.dll";
    pluginsdir.setNameFilters( namefilters );
#endif

#ifdef Q_WS_X11
    pluginsdir.setPath ( *it );
    QStringList namefilters;
    namefilters << "*.so";
    pluginsdir.setNameFilters ( namefilters );
#endif

#ifdef Q_OS_MAC
    pluginsdir.setPath ( *it );
// pluginsdir.setNameFilter("*.dylib");
#endif

    QStringList pluginslist=pluginsdir.entryList();

    QStringList::iterator mt;
    typedef  Plugin* ( *pf ) ();
    std::cout << "Searching on " << ( *it ).toStdString() << std::endl;
    //lib.setLoadHints()
    for ( mt=pluginslist.begin();mt!=pluginslist.end();mt++ )
    {
      QString filename= ( *it ) +QDir::separator() + ( *mt );
      pf function= ( pf ) QLibrary::resolve ( filename,"GetPlugin" );
      if ( function )
      {

        Plugin* pJ=function();
        m_private->m_plugins.push_back ( pJ );
        std::cout << pJ->Name() << " plugin loaded."  << std::endl;
        pluginfound=true;
      }

    }
  }
  if ( !pluginfound )
  {

    std::cout << "Any QryoMol plugin seems to be installed on your system" << std::endl;

  }

}



