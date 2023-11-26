/*****************************************************************************************
                            pluginloader.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QRYOMOLPLUGINLOADER_H
#define QRYOMOLPLUGINLOADER_H

#include <vector>

#include <QStringList>

#include "export.h"

namespace kryomol
{
class PluginLoaderPrivate;
class Plugin;

/**
  \brief A class to handle loading of plugins in QryoMol
  Thjis class manage the loading of qryomol Plugins in a platform independent way.
  By default application search for plugins inside the plugins directory of the distribution. Additional directories
  can be added using the QRYOMOL_PLUGINS_DIRS environment variable
  */
class KRYOMOL_API PluginLoader
{
public:
    PluginLoader();
    ~PluginLoader();
    std::vector<kryomol::Plugin*>& Plugins();
    const QStringList& Dirs() const;
    void Load();

private:
    void InitializePluginDirs();

private:
    PluginLoaderPrivate* m_private;

};

}

#endif
