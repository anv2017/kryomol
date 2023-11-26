/*****************************************************************************************
                            plugin.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "plugin.h"
#include "glvisor.h"


class kryomol::PluginPrivate
{
  public:
    PluginPrivate() : m_stream ( NULL ), m_world ( NULL )
    {}

    ~PluginPrivate()
    {}

    /*Pointer to the graphical elements*/
    QWidget* m_ui;
    World* m_world;
    QDockArea* m_dockarea;
    QMainWindow* m_mainwindow;
    /*associated stream to put output*/
    std::ostream* m_stream;
};

using namespace kryomol;



/**
 \brief Empty constructor
 Constructs an empty Plugin object
 */

Plugin:: Plugin()
{
  _d = new PluginPrivate();
}

/** 
\brief Virtual destructor
Reimplement in the derived classes if you need to do
further cleaning*/

Plugin::~Plugin()
{
  delete _d;
}

/** \brief Plugin initialization
This function is called for each plugin after runtime loading.inside the main application
Any initialization depending on molecular properties should be done here*/
void Plugin::Initialize()
{

}

/** \return graphical interface associated with the plugin.
In most cases the associated graphical interface will be directly derived from QWidget. Further graphical elements as additional QDockWidgets can be childs of this widget.
After runtime plugin loading the QryoMolMainWindow will store the widge as a child of a QDockWidget. Each widget plugin is hide or shown
using the Analysis menu tool*/
QWidget* Plugin::UI()
{
  return _d->m_ui;
}

/**Associate the simulation world w with the plugin*/
void Plugin::SetWorld ( World* w )
{
  _d->m_world=w;
}

/** Return a pointer to the associated QDockArea*/
QDockArea* Plugin::DockArea()
{
  return _d->m_dockarea;
}
/** Associate a QDockArea with the plugin
    In most cases this will be QDockArea of QryoMol main window*/
void Plugin::SetDockArea ( QDockArea* area )
{
  _d->m_dockarea=area;
}

/** \return a pointer to the associated simulation world*/
World* Plugin::GetWorld()
{
  return _d->m_world;
}

/** \return a const pointer to the associated simulation world*/
const World* Plugin::GetWorld() const
{
  return _d->m_world;
}

/** \return the name of the plugin. This name will be shown in the QryoMol Analysis menu and the
     plugins combobox*/
std::string Plugin::Name() const
{
  return "";
}

/** \return a plain text (ascii) brief description of the plugin*/
std::string Plugin::Description() const
{
  return "";
}

/** \return a detalied description of the plugin as ascii html formatted text*/
std::string Plugin::HTMLDescription() const
{
  return "";
}

/** \return a pointer to the associated main window*/
QMainWindow* Plugin::MainWindow() const
{
  return _d->m_mainwindow;
}

/** \brief associate the QryoMolMainWindow with the plugin
     This function is called for all plugins after runtime loading
     */
void Plugin::SetMainWindow ( QMainWindow* w )
{
  _d->m_mainwindow=w;
}

void Plugin::SetUI ( QWidget* w )
{
  _d->m_ui=w;
}

/** \brief attach a stream to the plugin for ouput of results
This function alllows to print the ouput of computational results inside
your plugin in a selected stream. For instance, if you want to show the results on the 
standard ouput
\code
MyPlugin::Initialize()
{
  AttachStream(&std::cout)
  
}
\endcode
*/
void Plugin::AttachStream ( std::ostream* s )
{
  _d->m_stream=s;
}

/** \brief dettach the associated stream*/
void Plugin::DettachStream()
{
  _d->m_stream=NULL;
}

/** \return a reference to the attached stl stream*/
std::ostream& Plugin::Stream()
{
  return * ( _d->m_stream );

}

/** \return a reference the attached stl stream*/
const std::ostream& Plugin::Stream() const
{
  return * ( _d->m_stream );
}

/** \brief handling of selection events
     When the user picks over ther atom \a atom the ProcessSelection function
     is called for the running plugin. Overriding this function in your plugin allows to control
     the selection of indivudual atoms
*/
void Plugin::ProcessSelection ( int atom )
{

}
/** Do further rendering*/
void Plugin::Render()
{
}




