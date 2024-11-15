/*****************************************************************************************
                            kryomolmainwindow.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef MOBDEVWINDOW_H
#define MOBDEVWINDOW_H

#include <QMainWindow>
#include <QStackedWidget>
#include <QString>
#include <QToolBar>
#include <QActionGroup>
#include <QProgressBar>
#include <QPair>
#include <qmolecularlistcontrol.h>
#include <qpdbcontrol.h>
#include <qmeasurewidget.h>
#include <qorbitalwidget.h>

#include "glvisor.h"
#include "kryovisor.h"
#include "renderorbitals.h"
#include "parser.h"

namespace kryomol
{
class World;
class KryoMolApplication;
class GLVisor;
class Parser;
class ParserFactory;
class KryoVisor;
class RenderOrbitals;
}

class OrcaDialog;


class KryoMolMainWindow : public QMainWindow
{

     Q_OBJECT

public:
    explicit KryoMolMainWindow(QWidget *parent = 0);
    ~KryoMolMainWindow();

public slots:
    void OnOpenFile(QString filename);

private:
    enum filetype { PLAINFILE=0,UVFOLDER,IRFOLDER };
    void UpdateRecentFiles(const QPair<QString,int>& f=QPair<QString,int>() );
    void OpenFile();
    void OpenUVFolder(QString foldername);
    void OpenIRFolder(QString foldername);
    void Init();
    void InitToolBars();
    void InitGaussian();
    void FinishGaussian();
    void OpenGenericFile();
    void OpenGaussianFile();
    void Parse();
    void ParseGaussian();
    bool isGaussianFile();
    bool isGaussianFile(QString file);

signals:
    void ParseringFinished();
    void showDensity(bool);

private slots:
    void OnLoadButton();
    void OnLoadUVFolderAction();
    void OnLoadIRFolderAction();
    void OnOpenRecentFile();
    void OnOpenConfButton();
    void OnOpenPDBButton();
    void OnOpenPluginButton();
    void OnOpenOrbitalButton(bool);
    void OnChangeJob();
    void OnChangeJobN();
    void OnChangedAnalysisMenu(QAction* );
    void OnChangePlugin(int);
    void OnStartAnimation();
    void OnStopAnimation();
    void OnPreviousFrame();
    void OnNextFrame();
    void OnFirstFrame();
    void OnLastFrame();
    void OnNoneSelection();
    void OnMeasureDistances();
    void OnMeasureAngles();
    void OnMeasureDihedrals();
    void OnRotateBond();
    void OnShowMeasures();
    void OnDrawDensity(bool);
    void OnUpdateMeasures(size_t);
    void OnUpdateDensities(size_t);
    void OnUpdateOrbitals(size_t);
    void OnExportVectorGraphics();
    void OnExportRasterGraphics();
    void OnExportGeom();
    void OnExportGeomCurrent();
    void OnExportACES();
    void OnExportCPMD();
    void OnExportGaussian(bool whithessian=false);
    void OnExportGaussianForces(bool whithessian=true);
    void OnSetGraphMode(QAction* );
    void OnCenterVisiblePart();
    void OnCenterWholeMolecule();
    void OnAbout();
    void OnRunOrcaWidget();
    void OnProtonateTrigonalCenterAction();
    void OnProtonateTrigonalCenter(std::vector<size_t> selatoms);
    void OnUpdateFrameForOrbitals(size_t f);

private:
    void InitWidgets(kryomol::JobType,bool hasdensity,bool hasorbitals);
protected slots:
    void OnGraphModeChanged( kryomol::GLVisor::graphmode );
    void SetBondOrders();

private:
    QStackedWidget* m_stackedwidget;
    QToolBar* m_maintoolbar;
    QList<QToolBar*> m_navigationtoolbar;
    QString m_file;
    kryomol::KryoMolApplication* m_app;
    std::vector<QAction*> m_recentfileactions;
    QActionGroup* m_gmodegroup;
    QMenu* m_analysismenu;
    QMolecularListControl* m_mdcontrol;
    QPDBControl* m_pdbcontrol;
    QStackedWidget* m_uistack;
    QMeasureWidget* m_measures;
    QOrbitalWidget* m_orbitals;
    int m_currentplugin;
    kryomol::Parser* m_qparser;
    std::vector<kryomol::GLVisor*> m_visors;
    std::vector<kryomol::World*> m_worlds;
    QSize m_size;
    QActionGroup* m_measureActions;
    QActionGroup* m_jobs;
    OrcaDialog* m_orcawidget;
    bool m_hasdensity;
    bool m_hasorbitals;
    bool m_hasalphabeta;
    QTabWidget* m_tabwidget;

};

#endif // MOBDEVWINDOW_H
