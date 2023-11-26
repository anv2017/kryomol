/*****************************************************************************************
                            plugin.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef PLUGIN_H
#define PLUGIN_H


#include <string>
#include <vector>
#include <iostream>


class QPixmap;
class QWidget;
class QDockArea;
class QMainWindow;

#include "export.h"

namespace kryomol
{
class World;
class Molecule;
class PluginPrivate;

/**
\brief The base class for qryomol plugins

Each kryomol plugin, a runtime loadable shared library, should inherit from the Plugin class
\code
class MyPlugin : public qryomol::Plugin()
{
}
\endcode

In order to be load at runtime an GetPlugin function block should be include in the implementation
returning a pointer to the plugin
\code
Plugin* GetPlugin()
{
  return new MyPlugin();
}
\endcode
*/
class  KRYOMOL_API Plugin
{
public:
    Plugin();
    virtual ~ Plugin();
    virtual void Initialize();
    QWidget* UI();
    void SetUI ( QWidget* w );
    void SetWorld ( World* w );
    QDockArea* DockArea();
    void SetDockArea ( QDockArea* area );
    QMainWindow* MainWindow() const;
    void SetMainWindow ( QMainWindow* w );
    World* GetWorld();
    const World* GetWorld() const;
    virtual std::string Name() const =0;
    virtual std::string Description() const;
    virtual std::string HTMLDescription() const;
    virtual void ProcessSelection ( int atom );
    virtual void Render();
    void AttachStream ( std::ostream* s );
    void DettachStream();
    std::ostream& Stream();
    const std::ostream& Stream() const;
private:
    PluginPrivate* _d;

};
};

/** 
    \return a pointer to a new qryomol::Plugin object
    This function should be always implemented in any new plugin
    \code
    Plugin* GetPlugin()
    {
      return new MyPlugin;
    }
*/
extern "C"
{
KRYOMOL_EXPORT kryomol::Plugin* GetPlugin();

}



#endif
