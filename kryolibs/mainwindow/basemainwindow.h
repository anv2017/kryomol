/*****************************************************************************************
                            basemainwindow.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef BASEMAINWINDOW_H
#define BASEMAINWINDOW_H

#include <vector>

#include <QDragEnterEvent>
#include <QDropEvent>
#include <QKeyEvent>
#include <QMainWindow>

#include "glvisor.h"

class QStackedWidget;
class QMimeData;
class QMolecularListControl;
class QComboBox;
class QAction;
class QPDBControl;
class QActionGroup;
class QOrbitalWidget;


namespace kryomol
{
class Parser;
class Plugin;
class World;
class PluginLoader;
class KryoMolApplication;
class Scripter;
class ParserFactory;
};
class QCWTComboBox;

/** \brief The KryoMol main window
*/

class BaseMainWindow : public QMainWindow
{
  Q_OBJECT
public:
  BaseMainWindow(QWidget* parent=0, Qt::WindowFlags f=0);
  BaseMainWindow(const char* file, QWidget* parent=0, Qt::WindowFlags f=0);
  ~BaseMainWindow();
  /** \return a pointer to the current plugin (loaded as a shared library*/
  kryomol::Plugin* CurrentPlugin();
 /** Load plugins for the experiments*/
  void LoadPlugins();
protected: // Protected methods
  /** No descriptions */
  virtual void keyPressEvent(QKeyEvent* e);
  /** Parse the molecular structure file with path s*/
  void OpenFile();
  //void OpenFileGaussian();
  /** update the QSettings for recent opened files*/
  void UpdateRecentFiles();
  bool isGaussianFile();
  bool existDensity();
  /**Parse the streams*/
  void Parse();
  //void ParseGaussian();
private:
  void InitControls();
private slots:
  /**Print screen*/
  void OnPrintButton();
  /** Change to the ith plugin*/
  void OnChangePlugin(int i);
  /** About stuff*/
  virtual void OnAbout();
  /** Go to home page*/
  void OnNMRDevOnline();
  /** Help*/
  void OnHelp();
  /**Paste Molecule*/
  void OnPaste();
  void DropHandler(const QMimeData* src);
  void OnChangedPluginsList(int i); 
  void OnChangedAnalysisMenu(QAction* );
  void OnExportVectorGraphics();
  void OnExportRasterGraphics();
  void OnPreviousFrame();
  void OnNextFrame(); 
  void OnLastFrame();
  void OnFirstFrame();
  void OnOpenRecentFile();
  void OnCenterVisiblePart();
  void OnCenterWholeMolecule();
  void OnExportGeom();
  void OnSetGraphMode(QAction* );
  void OnShowCommandLine(bool );
  void OnPluginTriggered(QAction* );

protected slots:
  void OnGraphModeChanged( kryomol::GLVisorBase::graphmode );
  void SetBondOrders();

protected:
  /** absolute path to the opened file*/
  QString m_file;
  kryomol::Parser* m_qparser;
  /** number of the current plugin*/
  int m_currentplugin;
  QToolBar* m_maintoolbar;
  QToolBar* m_navigationtoolbar;
  QCWTComboBox* m_pluginscombo;
  std::vector<QAction*> m_recentfileactions;
  QActionGroup* m_gmodegroup;
  QMenu* m_analysismenu;
  QMolecularListControl* m_mdcontrol;
  QPDBControl* m_pdbcontrol;
  QStackedWidget* m_uistack;
  kryomol::KryoMolApplication* m_app;
  QActionGroup* m_pluginactions;
  /** pointer to the simulation world*/
  kryomol::World* m_world;
  kryomol::Scripter* m_scripter;
};

#endif
