/*****************************************************************************************
                            qryomolmainwindow.cpp  -  description
                             -------------------
******************************************************************************************/

#include <sstream>
#include <stdlib.h>
#include <iostream>

#include "kryomolmainwindow.h"
#include "qryomolapp.h"
#include "world.h"
#include "atom.h"
#include "glvisor.h"
#include "parserfactory.h"
#include "stringtools.h"
#include "xyzwriter.h"
#include "aceswriter.h"
#include "cpmdwriter.h"
#include "kryovisor.h"
#include "parsers.h"
#include "orcadialog.h"

#include <QStackedWidget>
#include <QToolBar>
#include <QAction>
#include <QIcon>
#include <QString>
#include <QFileDialog>
#include <QMenu>
#include <QWhatsThis>
#include <QSplitter>
#include <QSettings>
#include <QMessageBox>
#include <QImageWriter>
#include <QClipboard>
#include <QTextStream>
#include <QStringList>
#include <QMenuBar>
#include <QProgressDialog>
#include <QStatusBar>
#include <QDockWidget>
#include <QTextEdit>

#include "qryomolinfo.h"
#include "qjobwidget.h"
#include "qjobspwidget.h"
#include "qjoboptwidget.h"
#include "qjobdynwidget.h"
#include "qjobfreqwidget.h"
#include "qjobuvwidget.h"
#include "qconvwidget.h"
#include "qfreqwidget.h"
#include "qirwidget.h"
#include "qjcdrawing.h"
#include "qmeasurewidget.h"
#include "openingthread.h"

#ifdef __MINGW32__
#include <filesystem>
#endif

using namespace kryomol;

const int maxrecentfiles = 5;
const int fifosize = 20;

KryoMolMainWindow::KryoMolMainWindow(QWidget *parent) : QMainWindow(),
    m_orcawidget(nullptr) , m_hasdensity(false), m_hasorbitals(false),
    m_hasalphabeta(false)
{
    Init();
}

KryoMolMainWindow::~KryoMolMainWindow()
{
}

void KryoMolMainWindow::InitToolBars()
{

    QMenu* filemenu = new QMenu(tr ( "File" ),this );
    QMenuBar* mainmenu = menuBar();

    mainmenu->addMenu(filemenu);

    //File menu
    QAction*   loadAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/fileopen.png") ),tr ( "&Open Molecule" ),this );
    loadAction->setShortcut(tr("Ctrl+O"));
    loadAction->setStatusTip(tr("Open an existing file"));
    connect ( loadAction,SIGNAL ( triggered() ),this,SLOT ( OnLoadButton() ) );
    filemenu->addAction ( loadAction );


    //Add openfolder action. This is by now intended to load a folder with several ECD computations
    QMenu* openfoldermenu = new QMenu("Open Folder");
    filemenu->addMenu(openfoldermenu);
    QAction*   loadUVFolderAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/fileopen.png") ),tr ( "&Open UV/ECD Folder " ),this );
    //loadFolderAction->setShortcut(tr("Ctrl+F"));
    loadUVFolderAction->setStatusTip(tr("Open many UV/ECD files from a folder"));
    connect ( loadUVFolderAction,SIGNAL ( triggered() ),this,SLOT ( OnLoadUVFolderAction() ) );
    openfoldermenu->addAction ( loadUVFolderAction );
    QAction*   loadIRFolderAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/fileopen.png") ),tr ( "&Open IR/VCD Folder " ),this );
    //loadFolderAction->setShortcut(tr("Ctrl+F"));
    loadIRFolderAction->setStatusTip(tr("Open many IR/VCD files from a folder"));
    connect ( loadIRFolderAction,SIGNAL ( triggered() ),this,SLOT ( OnLoadIRFolderAction() ) );
    openfoldermenu->addAction ( loadIRFolderAction );

    QAction* quitAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/exit.png") ),tr ( "&Quit" ),this );
    quitAction->setShortcut(tr("Ctrl+Q"));
    quitAction->setStatusTip(tr("Quit KryoMol"));
    connect ( quitAction,SIGNAL ( triggered() ),this,SLOT ( close() ) );
    filemenu->addAction ( quitAction );
    QMenu* recentFilesMenu = new QMenu ( tr ( "Recent Files" ),this );
    filemenu->addMenu ( recentFilesMenu );

    m_recentfileactions.resize ( maxrecentfiles );
    for ( int i = 0; i < maxrecentfiles; ++i )
    {
        m_recentfileactions.at ( i ) = new QAction ( this );
        m_recentfileactions.at ( i )->setVisible ( false );
        connect ( m_recentfileactions.at ( i ), SIGNAL ( triggered() ),this, SLOT ( OnOpenRecentFile() ) );
        recentFilesMenu->addAction ( m_recentfileactions.at ( i ) );
    }
    UpdateRecentFiles();

    //Edit Menu
    QMenu* editmenu = new QMenu ( tr ( "Edit" ),this );
    QMenu* exportgraphicsmenu = new QMenu ( tr ( "Export Graphics" ),this );
    editmenu->addMenu ( exportgraphicsmenu );
    QMenu* exportgeommenu = new QMenu( tr("Export Geometry"),this);
    editmenu->addMenu( exportgeommenu );
    QMenu* exportXYZmenu = new QMenu( tr("Export XYZ Coordinates"),this);
    exportgeommenu->addMenu( exportXYZmenu );

    //Help menu
    QMenu* helpmenu = new QMenu ( tr ( "Help" ),this );
    QAction* aboutAction = new QAction(tr("&About QryoMol"),this);
    connect ( aboutAction,SIGNAL ( triggered() ),this,SLOT ( OnAbout() ) );
    helpmenu->addAction ( aboutAction );
    mainmenu->addMenu(helpmenu);

}


void KryoMolMainWindow::Init()
{
    m_app=dynamic_cast<kryomol::KryoMolApplication*>qApp;;
    kryomol::BuildPeriodicTable();
    m_world = new kryomol::World(this);
    m_stackedwidget = new QStackedWidget (this);

    //Construct the visor and widgets

    m_world = new kryomol::World(this,World::glvisor);

    m_mdcontrol = new QMolecularListControl(this);
    m_pdbcontrol= new QPDBControl(this);
    m_uistack = new QStackedWidget(this);
    m_measures = new QMeasureWidget(this);
    m_orbitals = new QOrbitalWidget(this);

    connect( m_world->Visor(), SIGNAL ( distance ( QString& ) ), m_measures, SLOT ( OnWriteDistance ( QString& ) ) );
    connect( m_world->Visor(), SIGNAL ( angle ( QString& ) ), m_measures, SLOT ( OnWriteAngle ( QString& ) ) );
    connect( m_world->Visor(), SIGNAL ( dihedral ( QString& ) ), m_measures, SLOT ( OnWriteDihedral ( QString& ) ) );
    connect(m_world, SIGNAL ( currentFrame ( size_t ) ), this, SLOT ( OnUpdateMeasures ( size_t )) );
    connect(m_measures,SIGNAL(clearAll()),m_world->Visor(), SLOT(OnClearMeasures()));
    connect(m_measures,SIGNAL(distanceChange(int)),m_world->Visor(),SLOT(OnDistanceChange(int)));
    connect(m_measures,SIGNAL(angleChange(int)),m_world->Visor(),SLOT(OnAngleChange(int)));
    connect(m_measures,SIGNAL(dihedralChange(int)),m_world->Visor(),SLOT(OnDihedralChange(int)));
    connect(m_measures, SIGNAL(showDistances(bool)),m_world->Visor(),SLOT(OnShowDistances(bool)));

    m_size = this->size();

    QSplitter* q = new QSplitter(this);
    q->addWidget(m_world->Visor());
    q->addWidget(m_pdbcontrol);
    q->addWidget(m_mdcontrol);
    q->addWidget(m_uistack);
    q->addWidget(m_measures);
    q->addWidget(m_orbitals);

    QList<int> list = q->sizes();
    list[0] = m_size.width();
    list[1] = 0;
    list[2] = 0;
    list[3] = 0;
    list[4] = 0;
    list[5] = 0;
    q->setSizes(list);

    //Initially only the visor is showed
    m_pdbcontrol->hide();
    m_mdcontrol->hide();
    m_uistack->hide();
    m_measures->hide();
    m_orbitals->hide();

    m_stackedwidget->addWidget(q);

    setCentralWidget(m_stackedwidget);

    //Construct the tool bars

    //----------------- Main Tool Bar -------------------

    //Actions
    QAction* openConfAction = new QAction("Conf",this);
    connect(openConfAction,SIGNAL(triggered()),this,SLOT(OnOpenConfButton()));

    QAction* openConfActionMenu = new QAction("Show Conformers",this);
    openConfActionMenu->setCheckable ( true );
    openConfActionMenu->setChecked ( false );
    connect(openConfActionMenu,SIGNAL(triggered()),this,SLOT(OnOpenConfButton()));

    QAction* openPDBAction = new QAction("PDB",this);
    connect(openPDBAction,SIGNAL(triggered()),this,SLOT(OnOpenPDBButton()));

    QAction* openPDBActionMenu = new QAction("Show PDB",this);
    openPDBActionMenu->setCheckable ( true );
    openPDBActionMenu->setChecked ( false );
    connect(openPDBActionMenu,SIGNAL(triggered()),this,SLOT(OnOpenPDBButton()));

    QAction* openPluginAction = new QAction("Plugins",this);
    connect(openPluginAction,SIGNAL(triggered()),this,SLOT(OnOpenPluginButton()));

    QAction* openPluginActionMenu = new QAction("Show Plugins",this);
    openPluginActionMenu->setCheckable ( true );
    openPluginActionMenu->setChecked ( false );
    connect(openPluginActionMenu,SIGNAL(triggered()),this,SLOT(OnOpenPluginButton()));

    //QPushButton* showMeasures = new QPushButton(tr("    Show \n Measurements"),this);
    //connect(showMeasures,SIGNAL(clicked()), this, SLOT(OnShowMeasures()));

    QAction* showMeasuresAction = new QAction(tr("Show Measurements"),this);
    showMeasuresAction->setCheckable ( true );
    showMeasuresAction->setChecked ( false );
    connect(showMeasuresAction,SIGNAL(triggered()), this, SLOT(OnShowMeasures()));


    //File menu
    QMenu* filemenu = new QMenu ( tr ( "&File" ),this );
    QAction*   loadAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/fileopen.png") ),tr ( "&Open Molecule" ),this );
    loadAction->setShortcut(tr("Ctrl+O"));
    loadAction->setStatusTip(tr("Open an existing file"));
    connect ( loadAction,SIGNAL ( triggered() ),this,SLOT ( OnLoadButton() ) );
    filemenu->addAction ( loadAction );

    QAction* quitAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/exit.png") ),tr ( "&Quit" ),this );
    quitAction->setShortcut(tr("Ctrl+Q"));
    quitAction->setStatusTip(tr("Quit KryoMol"));
    connect ( quitAction,SIGNAL ( triggered() ),this,SLOT ( close() ) );
    filemenu->addAction ( quitAction );
    QMenu* recentFilesMenu = new QMenu ( tr ( "Recent Files" ),this );
    filemenu->addMenu ( recentFilesMenu );

    m_recentfileactions.resize ( maxrecentfiles );
    for ( int i = 0; i < maxrecentfiles; ++i )
    {
        m_recentfileactions.at ( i ) = new QAction ( this );
        m_recentfileactions.at ( i )->setVisible ( false );
        connect ( m_recentfileactions.at ( i ), SIGNAL ( triggered() ),this, SLOT ( OnOpenRecentFile() ) );
        recentFilesMenu->addAction ( m_recentfileactions.at ( i ) );
    }
    UpdateRecentFiles();


    //Edit Menu
    QMenu* editmenu = new QMenu ( tr ( "Edit" ),this );
    QMenu* exportgraphicsmenu = new QMenu ( tr ( "Export Graphics" ),this );
    editmenu->addMenu ( exportgraphicsmenu );
    QMenu* exportgeommenu = new QMenu( tr("Export Geometry"),this);
    editmenu->addMenu( exportgeommenu );
    QMenu* exportXYZmenu = new QMenu( tr("Export XYZ Coordinates"),this);
    exportgeommenu->addMenu( exportXYZmenu );

    QAction* enantiomerAction = new QAction( tr("Enantiomerize"),this);
    connect ( enantiomerAction,SIGNAL ( triggered() ),m_world->Visor(), SLOT ( OnEnantiomerize() ) );
    editmenu->addAction( enantiomerAction );

    QAction* exportRasterAction = new QAction ( tr ( "Export Raster Graphics" ),this );
    connect ( exportRasterAction,SIGNAL ( triggered() ),this, SLOT ( OnExportRasterGraphics() ) );
    exportgraphicsmenu->addAction ( exportRasterAction );

    QAction* exportVectorAction = new QAction ( tr ( "Export Vector Graphics" ),this );
    connect ( exportVectorAction,SIGNAL ( triggered() ),this, SLOT ( OnExportVectorGraphics() ) );
    exportgraphicsmenu->addAction ( exportVectorAction );

    QAction* exportACESAction = new QAction ( tr ( "Export CFOUR Input" ),this );
    connect ( exportACESAction,SIGNAL ( triggered() ),this, SLOT ( OnExportACES() ) );
    exportgeommenu->addAction ( exportACESAction );

    QAction* exportCPMDAction = new QAction ( tr ( "Export CPMD Input" ),this );
    connect ( exportCPMDAction,SIGNAL ( triggered() ),this, SLOT ( OnExportCPMD() ) );
    exportgeommenu->addAction ( exportCPMDAction );

    QAction* exportGaussianAction = new QAction ( tr ( "Export Gaussian Input" ),this );
    connect ( exportGaussianAction,SIGNAL ( triggered()),this, SLOT ( OnExportGaussian() ) );
    exportgeommenu->addAction ( exportGaussianAction );

    QAction* exportXYZAction = new QAction( tr("All"),this);
    connect(exportXYZAction,SIGNAL( triggered() ),this,SLOT(OnExportGeom()));
    exportXYZmenu->addAction(exportXYZAction);

    QAction* exportXYZcurrentAction = new QAction( tr("Current"),this);
    connect(exportXYZcurrentAction,SIGNAL( triggered() ),this,SLOT(OnExportGeomCurrent()));
    exportXYZmenu->addAction(exportXYZcurrentAction);

    QMenu* extratoolsmenu= new QMenu(tr ( "Extra tools" ),this );
    editmenu->addMenu(extratoolsmenu);
    QAction* protonateTrigonalCenterAction = new QAction( tr("Protonate trigonal center"),this);
    connect(protonateTrigonalCenterAction,SIGNAL(triggered()),this,SLOT(OnProtonateTrigonalCenterAction()));
    QAction* orcaJobAction = new QAction( tr("Run Orca job") );
    connect(orcaJobAction, SIGNAL(triggered()),this,SLOT(OnRunOrcaWidget()));
    extratoolsmenu->addAction(protonateTrigonalCenterAction);
    extratoolsmenu->addAction(orcaJobAction);



    //View menu
    QMenu* viewmenu= new QMenu( tr("View"),this);
    viewmenu->setTearOffEnabled(true);
    QMenu* cameramenu= new QMenu( tr("Center"),this);
    viewmenu->addMenu(cameramenu);

    QAction* shownumbers= new QAction ( tr("Show atom numbers"),this );
    shownumbers->setCheckable ( true );
    shownumbers->setChecked ( m_world->Visor()->ShowNumbers() );
    connect ( shownumbers,SIGNAL ( toggled ( bool ) ),m_world->Visor(),SLOT ( OnShowNumbers ( bool ) ) );
    viewmenu->addAction ( shownumbers );

    QAction* showsymbols= new QAction ( tr("Show atom symbols") ,this );
    showsymbols->setCheckable ( true );
    showsymbols->setChecked ( m_world->Visor()->ShowSymbols() );
    connect ( showsymbols,SIGNAL ( toggled ( bool ) ),m_world->Visor(), SLOT ( OnShowSymbols ( bool ) ) );
    viewmenu->addAction ( showsymbols );

    QAction* showpdbinfo= new QAction ( tr("Show PDB Info") ,this );
    showpdbinfo->setCheckable ( true );
    showpdbinfo->setChecked ( m_world->Visor()->ShowPDBInfo() );
    connect ( showpdbinfo,SIGNAL ( toggled ( bool ) ),m_world->Visor(), SLOT ( OnShowPDBInfo ( bool ) ) );
    viewmenu->addAction ( showpdbinfo );

    QAction* showdipole= new QAction ( tr("Show Dipole Moment") ,this );
    showdipole->setCheckable ( true );
    showdipole->setChecked ( m_world->Visor()->ShowDipole() );
    connect ( showdipole,SIGNAL ( toggled ( bool ) ),m_world->Visor(), SLOT ( OnShowDipole ( bool ) ) );
    viewmenu->addAction ( showdipole );

    QAction* showcell= new QAction ( tr("Show Cell") ,this );
    showcell->setCheckable ( true );
    showcell->setChecked ( m_world->Visor()->ShowCell() );
    connect ( showcell,SIGNAL ( toggled ( bool ) ),m_world->Visor(), SLOT ( OnShowCell ( bool ) ) );
    viewmenu->addAction ( showcell );

    QAction* showdensity= new QAction ( tr("Show Electronic Density") ,this );
    if ( m_hasdensity )
    {
        showdensity->setEnabled(true);
        showdensity->setCheckable ( true );
        showdensity->setChecked ( m_world->Visor()->ShowDensity() );
        connect ( showdensity,SIGNAL ( toggled ( bool ) ),this, SLOT ( OnOpenOrbitalButton ( bool ) ) );
        connect ( m_orbitals, SIGNAL ( showDensity ( bool ) ), this, SLOT ( OnShowDensity ( bool ) ) );
        connect ( this, SIGNAL ( drawDensity ( bool ) ), m_world->Visor(), SLOT ( OnDrawDensity ( bool ) ) );
    }
    else
        showdensity->setEnabled(false);

    viewmenu->addAction ( showdensity );

    QAction* showaxis= new QAction ( tr("Show Coordinate Axis") ,this );
    showaxis->setCheckable ( true );
    showaxis->setChecked ( m_world->Visor()->ShowAxis() );
    connect ( showaxis,SIGNAL ( toggled ( bool ) ),m_world->Visor(), SLOT ( OnShowAxis ( bool ) ) );
    viewmenu->addAction ( showaxis );

    QAction* changebackground = new QAction ( "Change background color",this );
    connect ( changebackground,SIGNAL ( triggered() ),m_world->Visor(),SLOT ( OnChangeBackground() ) );
    viewmenu->addAction ( changebackground );
    QMenu* rendermenu= viewmenu->addMenu("Rendering Mode");
    QAction* wfaction=rendermenu->addAction ( "Wireframe" );
    wfaction->setData(QVariant(0));
    wfaction->setCheckable(true);
    QAction* sticksaction=rendermenu->addAction ( tr("Sticks") );
    sticksaction->setData(QVariant(1));
    sticksaction->setCheckable(true);
    QAction* cpkaction=rendermenu->addAction ( "CPK" );
    cpkaction->setData(QVariant(2));
    cpkaction->setCheckable(true);
    m_gmodegroup= new QActionGroup(this);
    m_gmodegroup->addAction(wfaction);
    m_gmodegroup->addAction(sticksaction);
    m_gmodegroup->addAction(cpkaction);
    m_gmodegroup->setExclusive(true);
    wfaction->setChecked(true);

    connect ( rendermenu,SIGNAL ( triggered ( QAction* ) ),this,SLOT ( OnSetGraphMode ( QAction* ) ) );

    QAction* cameraOnWholeAction = new QAction( tr ("Whole molecule"),this );
    connect (cameraOnWholeAction,SIGNAL( triggered() ),this,SLOT(OnCenterWholeMolecule()));
    QAction* cameraOnVisibleAction = new QAction( tr("Visible part"),this);
    connect( cameraOnVisibleAction,SIGNAL( triggered() ),this,SLOT(OnCenterVisiblePart()));
    cameramenu->addAction(cameraOnWholeAction);
    cameramenu->addAction(cameraOnVisibleAction);

    //Window menu
    QMenu* windowmenu = new QMenu(tr ("Window"),this );
    windowmenu->addAction(openConfActionMenu);
    windowmenu->addAction(openPDBActionMenu);
    windowmenu->addAction(openPluginActionMenu);
    windowmenu->addAction(showMeasuresAction);


    QMenuBar* mainmenu = menuBar();
    //Menus
    mainmenu->addMenu(filemenu);
    mainmenu->addMenu (editmenu);
    mainmenu->addMenu(viewmenu);
    mainmenu->addMenu(windowmenu);


    //---------------- Navigation Tool Bar ---------------

    QToolBar* navigationtoolbar = new QToolBar(this);

    //Actions for frames
    QAction* firstFrameAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/start.png") ),tr ( "First Frame" ),this );
    connect ( firstFrameAction,SIGNAL ( triggered() ),this,SLOT ( OnFirstFrame() ) );

    QAction* lastFrameAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/finish.png" ) ), tr ( "Last Frame" ),this );
    connect ( lastFrameAction,SIGNAL ( triggered() ),this,SLOT ( OnLastFrame() ) );

    QAction* nextFrameAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/forward.png") ),tr ( "Next Frame" ),this );
    connect ( nextFrameAction,SIGNAL ( triggered() ),this,SLOT ( OnNextFrame() ) );

    QAction* previousFrameAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/back.png") ),tr ( "Previous Frame" ),this );
    connect ( previousFrameAction,SIGNAL ( triggered() ),this,SLOT ( OnPreviousFrame() ) );

    navigationtoolbar->addAction ( firstFrameAction );
    navigationtoolbar->addAction ( lastFrameAction );
    navigationtoolbar->addAction ( previousFrameAction );
    navigationtoolbar->addAction ( nextFrameAction );

    //Actions for measures -> QActionGroup
    m_measureActions = new QActionGroup ( this );
    QAction* noneAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/select.png") ),tr ( "Plain Selection" ),this );
    noneAction->setActionGroup ( m_measureActions );
    connect ( noneAction,SIGNAL ( triggered() ),this,SLOT ( OnNoneSelection() ) );
    QAction* distanceAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/distance.png") ),tr ( "Measure Distances" ),this );
    distanceAction->setActionGroup ( m_measureActions );
    connect ( distanceAction,SIGNAL ( triggered() ),this,SLOT ( OnMeasureDistances () ) );

    QAction* angleAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/angle.png") ), tr ( "Measure Angles" ),this );
    connect ( angleAction,SIGNAL ( triggered() ),this,SLOT ( OnMeasureAngles () ) );
    angleAction->setActionGroup ( m_measureActions );
    QAction* dihedralAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/dihedral.png") ),tr ( "Measure Dihedrals" ),this );
    dihedralAction->setActionGroup ( m_measureActions );
    connect ( dihedralAction,SIGNAL ( triggered() ),this,SLOT ( OnMeasureDihedrals () ) );

    QAction* rotateBondAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/dihedral.png") ),tr ( "Rotate Bond" ),this );
    rotateBondAction->setActionGroup ( m_measureActions );
    connect ( rotateBondAction,SIGNAL ( triggered() ),this,SLOT ( OnRotateBond () ) );

    noneAction->setCheckable ( true );
    distanceAction->setCheckable ( true );
    angleAction->setCheckable ( true );
    dihedralAction->setCheckable ( true );
    rotateBondAction->setCheckable(true);
    m_measureActions->setExclusive ( true );

    navigationtoolbar->addAction ( noneAction );
    navigationtoolbar->addAction ( distanceAction );
    navigationtoolbar->addAction ( angleAction );
    navigationtoolbar->addAction ( dihedralAction );
    navigationtoolbar->addAction ( rotateBondAction );
    noneAction->setChecked ( true );

    m_navigationtoolbar.push_front(navigationtoolbar);

    //Add the navigation tool bar to the main window
    for(int i=0;i<m_navigationtoolbar.size();++i)
    {
        addToolBar(Qt::RightToolBarArea, m_navigationtoolbar.at(i));
        m_navigationtoolbar.at(i)->show();
    }

    update();

}


void KryoMolMainWindow::InitGaussian()
{
    //Elimate all widgets defined before
    delete m_stackedwidget;
    m_stackedwidget = new QStackedWidget (this);
    delete m_orcawidget;
    m_orcawidget=nullptr;

    m_size = this->size();

}

void KryoMolMainWindow::FinishGaussian()
{
    //elimate all tool bars defined before
    menuBar()->clear();


    for (int i=0;i<m_navigationtoolbar.size();++i)
    {
        removeToolBar(m_navigationtoolbar.at(i));
    }
    m_navigationtoolbar.clear();

    //Construct the tool bars

    //----------------- Main Tool Bar -------------------

    //Actions
    //QPushButton* changeJob = new QPushButton (tr("Change\n    Job"), this);
    //connect (changeJob, SIGNAL (clicked()),this,SLOT (OnChangeJob()));

    //QPushButton* showMeasures = new QPushButton(tr("    Show \n Measurements"),this);
    //connect(showMeasures,SIGNAL(clicked()), this, SLOT(OnShowMeasures()));

    QAction* showMeasuresAction = new QAction(tr("Show Measurements"),this);
    showMeasuresAction->setCheckable(true);
    showMeasuresAction->setChecked(false);
    connect(showMeasuresAction,SIGNAL ( triggered() ), this, SLOT(OnShowMeasures()));


    //File menu
    QMenu* filemenu = new QMenu ( tr ( "&File" ),this );
    QAction*   loadAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/fileopen.png") ),tr ( "&Open Molecule" ),this );
    loadAction->setShortcut(tr("Ctrl+O"));
    loadAction->setStatusTip(tr("Open an existing file"));
    connect ( loadAction,SIGNAL ( triggered() ),this,SLOT ( OnLoadButton() ) );
    filemenu->addAction ( loadAction );

    QAction* quitAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/exit.png") ),tr ( "&Quit" ),this );
    quitAction->setShortcut(tr("Ctrl+Q"));
    quitAction->setStatusTip(tr("Quit Qryomol"));
    connect ( quitAction,SIGNAL ( triggered() ),this,SLOT ( close() ) );
    filemenu->addAction ( quitAction );
    QMenu* recentFilesMenu = new QMenu ( tr ( "Recent Files" ),this );
    filemenu->addMenu ( recentFilesMenu );

    m_recentfileactions.resize ( maxrecentfiles );
    for ( int i = 0; i < maxrecentfiles; ++i )
    {
        m_recentfileactions.at ( i ) = new QAction ( this );
        m_recentfileactions.at ( i )->setVisible ( false );
        connect ( m_recentfileactions.at ( i ), SIGNAL ( triggered() ),this, SLOT ( OnOpenRecentFile() ) );
        recentFilesMenu->addAction ( m_recentfileactions.at ( i ) );
    }
    UpdateRecentFiles();


    //Edit Menu
    QMenu* editmenu = new QMenu ( tr ( "Edit" ),this );
    QMenu* exportgraphicsmenu = new QMenu ( tr ( "Export Graphics" ),this );
    editmenu->addMenu ( exportgraphicsmenu );
    QMenu* exportgeommenu = new QMenu( tr("Export Geometry"),this);
    editmenu->addMenu( exportgeommenu );
    QMenu* exportXYZmenu = new QMenu( tr("Export XYZ Coordinates"),this);
    exportgeommenu->addMenu( exportXYZmenu );
    QMenu* exportGaussianmenu = new QMenu( tr("Export Gaussian Input"),this);
    exportgeommenu->addMenu( exportGaussianmenu );

    QAction* enantiomerAction = new QAction( tr("Enantiomerize"),this);
    connect ( enantiomerAction,SIGNAL ( triggered() ),m_world->Visor(), SLOT ( OnEnantiomerize() ) );
    editmenu->addAction( enantiomerAction );

    QAction* exportRasterAction = new QAction ( tr ( "Export Raster Graphics" ),this );
    connect ( exportRasterAction,SIGNAL ( triggered() ),this, SLOT ( OnExportRasterGraphics() ) );
    exportgraphicsmenu->addAction ( exportRasterAction );

    QAction* exportVectorAction = new QAction ( tr ( "Export Vector Graphics" ),this );
    connect ( exportVectorAction,SIGNAL ( triggered() ),this, SLOT ( OnExportVectorGraphics() ) );
    exportgraphicsmenu->addAction ( exportVectorAction );

    QAction* exportACESAction = new QAction ( tr ( "Export CFOUR Input" ),this );
    connect ( exportACESAction,SIGNAL ( triggered() ),this, SLOT ( OnExportACES() ) );
    exportgeommenu->addAction ( exportACESAction );

    QAction* exportCPMDAction = new QAction ( tr ( "Export CPMD Input" ),this );
    connect ( exportCPMDAction,SIGNAL ( triggered() ),this, SLOT ( OnExportCPMD() ) );
    exportgeommenu->addAction ( exportCPMDAction );

    QAction* exportGaussianAction = new QAction ( tr ( "Geometry" ),this );
    connect ( exportGaussianAction,SIGNAL ( triggered() ),this, SLOT ( OnExportGaussian() ) );
    exportGaussianmenu->addAction ( exportGaussianAction );

    QJobWidget* current = ( static_cast<QJobWidget*> ( m_stackedwidget->currentWidget()));
    if ( current )
    {
        if (current->findChild<KryoVisorFreq*>())
        {
            QAction* exportGaussianForcesAction = new QAction ( tr ( "Geometry && Forces"),this );
            connect ( exportGaussianForcesAction,SIGNAL ( triggered() ),this, SLOT ( OnExportGaussianForces() ) );
            exportGaussianmenu->addAction ( exportGaussianForcesAction );
        }
    }

    QAction* exportNICSGridGaussianAction = new QAction ( tr ( "Geometry and NICS Grid" ),this );
    connect ( exportGaussianAction,SIGNAL ( triggered() ),this, SLOT ( OnExportGaussianNICSGrid() ) );
    exportGaussianmenu->addAction ( exportNICSGridGaussianAction );

    QAction* exportXYZAction = new QAction( tr("All"),this);
    connect(exportXYZAction,SIGNAL( triggered() ),this,SLOT(OnExportGeom()));
    exportXYZmenu->addAction(exportXYZAction);

    QAction* exportXYZcurrentAction = new QAction( tr("Current"),this);
    connect(exportXYZcurrentAction,SIGNAL( triggered() ),this,SLOT(OnExportGeomCurrent()));
    exportXYZmenu->addAction(exportXYZcurrentAction);

    QMenu* extratoolsmenu= new QMenu(tr ( "Extra tools" ),this );
    editmenu->addMenu(extratoolsmenu);
    QAction* protonateTrigonalCenterAction = new QAction( tr("Protonate trigonal center"),this);
    connect(protonateTrigonalCenterAction,SIGNAL(triggered()),this,SLOT(OnProtonateTrigonalCenterAction()));
    QAction* orcaJobAction = new QAction( tr("Run Orca job") );
    connect(orcaJobAction, SIGNAL(triggered()),this,SLOT(OnRunOrcaWidget()));
    extratoolsmenu->addAction(protonateTrigonalCenterAction);
    extratoolsmenu->addAction(orcaJobAction);




    //View menu
    QMenu* viewmenu= new QMenu( tr("View"),this);
    viewmenu->setTearOffEnabled(true);
    QMenu* cameramenu= new QMenu( tr("Center"),this);
    viewmenu->addMenu(cameramenu);

    QAction* shownumbers= new QAction ( tr("Show atom numbers"),this );
    shownumbers->setCheckable ( true );
    shownumbers->setChecked ( m_world->Visor()->ShowNumbers() );
    connect ( shownumbers,SIGNAL ( toggled ( bool ) ),m_world->Visor(),SLOT ( OnShowNumbers ( bool ) ) );
    viewmenu->addAction ( shownumbers );

    QAction* showsymbols= new QAction ( tr("Show atom symbols") ,this );
    showsymbols->setCheckable ( true );
    showsymbols->setChecked ( m_world->Visor()->ShowSymbols() );
    connect ( showsymbols,SIGNAL ( toggled ( bool ) ),m_world->Visor(), SLOT ( OnShowSymbols ( bool ) ) );
    viewmenu->addAction ( showsymbols );

    QAction* showpdbinfo= new QAction ( tr("Show PDB Info") ,this );
    showpdbinfo->setCheckable ( true );
    showpdbinfo->setChecked ( m_world->Visor()->ShowPDBInfo() );
    connect ( showpdbinfo,SIGNAL ( toggled ( bool ) ),m_world->Visor(), SLOT ( OnShowPDBInfo ( bool ) ) );
    viewmenu->addAction ( showpdbinfo );

    QAction* showdipole= new QAction ( tr("Show Dipole Moment") ,this );
    showdipole->setCheckable ( true );
    showdipole->setChecked ( m_world->Visor()->ShowDipole() );
    connect ( showdipole,SIGNAL ( toggled ( bool ) ),m_world->Visor(), SLOT ( OnShowDipole ( bool ) ) );
    viewmenu->addAction ( showdipole );

    QAction* showcell= new QAction ( tr("Show Cell") ,this );
    showcell->setCheckable ( true );
    showcell->setChecked ( m_world->Visor()->ShowCell() );
    connect ( showcell,SIGNAL ( toggled ( bool ) ),m_world->Visor(), SLOT ( OnShowCell ( bool ) ) );
    viewmenu->addAction ( showcell );

    QAction* showdensity= new QAction ( tr("Show Electronic Density") ,this );
    if ( m_hasdensity )
    {
        showdensity->setEnabled(true);
        showdensity->setCheckable ( true );
        showdensity->setChecked ( m_world->Visor()->ShowDensity() );
        connect ( showdensity,SIGNAL ( toggled ( bool ) ),this, SLOT ( OnOpenOrbitalButton ( bool ) ) );
    }
    else
        showdensity->setEnabled(false);
    viewmenu->addAction ( showdensity );

    QAction* showaxis= new QAction ( tr("Show Coordinate Axis") ,this );
    showaxis->setCheckable ( true );
    showaxis->setChecked ( m_world->Visor()->ShowAxis() );
    connect ( showaxis,SIGNAL ( toggled ( bool ) ),m_world->Visor(), SLOT ( OnShowAxis ( bool ) ) );
    viewmenu->addAction ( showaxis );

    QAction* changebackground = new QAction ( "Change background color",this );
    connect ( changebackground,SIGNAL ( triggered() ),m_world->Visor(),SLOT ( OnChangeBackground() ) );
    viewmenu->addAction ( changebackground );
    QMenu* rendermenu= viewmenu->addMenu("Rendering Mode");
    QAction* wfaction=rendermenu->addAction ( "Wireframe" );
    wfaction->setData(QVariant(0));
    wfaction->setCheckable(true);
    QAction* sticksaction=rendermenu->addAction ( tr("Sticks") );
    sticksaction->setData(QVariant(1));
    sticksaction->setCheckable(true);
    QAction* cpkaction=rendermenu->addAction ( "CPK" );
    cpkaction->setData(QVariant(2));
    cpkaction->setCheckable(true);
    m_gmodegroup= new QActionGroup(this);
    m_gmodegroup->addAction(wfaction);
    m_gmodegroup->addAction(sticksaction);
    m_gmodegroup->addAction(cpkaction);
    m_gmodegroup->setExclusive(true);
    wfaction->setChecked(true);

    connect ( rendermenu,SIGNAL ( triggered ( QAction* ) ),this,SLOT ( OnSetGraphMode ( QAction* ) ) );

    QAction* cameraOnWholeAction = new QAction( tr ("Whole molecule"),this );
    connect (cameraOnWholeAction,SIGNAL( triggered() ),this,SLOT(OnCenterWholeMolecule()));
    QAction* cameraOnVisibleAction = new QAction( tr("Visible part"),this);
    connect( cameraOnVisibleAction,SIGNAL( triggered() ),this,SLOT(OnCenterVisiblePart()));
    cameramenu->addAction(cameraOnWholeAction);
    cameramenu->addAction(cameraOnVisibleAction);


    //Job menu
    QMenu* jobmenu = new QMenu(tr("Job"),this);
    m_jobs = new QActionGroup(this);
    for (int i=1; i<=m_stackedwidget->count();i++)
    {
        QString name,num;
        num.sprintf("%d",i);
        name = "Job "+num;
        QAction* jobAction = new QAction(name,this);
        jobAction->setCheckable(true);
        jobAction->setActionGroup(m_jobs);
        if (i == m_stackedwidget->currentIndex()+1)
            jobAction->setChecked(true);
        else
            jobAction->setChecked(false);
        connect(jobAction,SIGNAL ( triggered() ), this, SLOT(OnChangeJobN()));
        jobmenu->addAction(jobAction);
    }
    m_jobs->setExclusive(true);


    //Window menu
    QMenu* windowmenu = new QMenu  (tr("Window"),this);
    windowmenu->addAction(showMeasuresAction);

    //Menu Bar, mainmenu, if we use Mac or Windows <-> Vertical Tool Bar, m_maintoolbar, on the left if we use a mobile device

    QMenuBar* mainmenu = menuBar();
    //Menus
    mainmenu->addMenu(filemenu);
    mainmenu->addMenu (editmenu);
    mainmenu->addMenu(viewmenu);
    mainmenu->addMenu(jobmenu);
    mainmenu->addMenu(windowmenu);

    //---------------- Navigation Tool Bar ---------------

    QToolBar* navigationtoolbar = new QToolBar(this);

    //Actions for frames
    QAction* firstFrameAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/start.png") ),tr ( "First Frame" ),this );
    connect ( firstFrameAction,SIGNAL ( triggered() ),this,SLOT ( OnFirstFrame() ) );

    QAction* lastFrameAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/finish.png" ) ), tr ( "Last Frame" ),this );
    connect ( lastFrameAction,SIGNAL ( triggered() ),this,SLOT ( OnLastFrame() ) );

    QAction* nextFrameAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/forward.png") ),tr ( "Next Frame" ),this );
    connect ( nextFrameAction,SIGNAL ( triggered() ),this,SLOT ( OnNextFrame() ) );

    QAction* previousFrameAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/back.png") ),tr ( "Previous Frame" ),this );
    connect ( previousFrameAction,SIGNAL ( triggered() ),this,SLOT ( OnPreviousFrame() ) );

    QAction* startAnimationAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/player_play.png") ),tr ( "Start Animation" ),this );
    connect ( startAnimationAction,SIGNAL ( triggered() ),this,SLOT ( OnStartAnimation() ) );

    QAction* stopAnimationAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/player_stop.png") ),tr ( "Stop Animation" ),this );
    connect ( stopAnimationAction,SIGNAL ( triggered() ),this,SLOT ( OnStopAnimation() ) );

    navigationtoolbar->addAction ( firstFrameAction );
    navigationtoolbar->addAction ( lastFrameAction );
    navigationtoolbar->addAction ( previousFrameAction );
    navigationtoolbar->addAction ( nextFrameAction );
    navigationtoolbar->addAction ( startAnimationAction );
    navigationtoolbar->addAction ( stopAnimationAction );

    //Actions for measures -> QActionGroup
    m_measureActions = new QActionGroup ( this );
    QAction* noneAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/select.png") ),tr ( "Plain Selection" ),this );
    noneAction->setActionGroup ( m_measureActions );
    connect ( noneAction,SIGNAL ( triggered() ),this,SLOT ( OnNoneSelection() ) );
    QAction* distanceAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/distance.png") ),tr ( "Measure Distances" ),this );
    distanceAction->setActionGroup ( m_measureActions );
    connect ( distanceAction,SIGNAL ( triggered() ),this,SLOT ( OnMeasureDistances () ) );

    QAction* angleAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/angle.png") ), tr ( "Measure Angles" ),this );
    connect ( angleAction,SIGNAL ( triggered() ),this,SLOT ( OnMeasureAngles () ) );
    angleAction->setActionGroup ( m_measureActions );
    QAction* dihedralAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/dihedral.png") ),tr ( "Measure Dihedrals" ),this );
    dihedralAction->setActionGroup ( m_measureActions );
    connect ( dihedralAction,SIGNAL ( triggered() ),this,SLOT ( OnMeasureDihedrals () ) );
    QAction* rotateBondAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/rotatebond.png") ),tr ( "Rotate Bond" ),this );
    rotateBondAction->setActionGroup ( m_measureActions );
    connect ( rotateBondAction,SIGNAL ( triggered() ),this,SLOT ( OnRotateBond() ) );

    noneAction->setCheckable ( true );
    distanceAction->setCheckable ( true );
    angleAction->setCheckable ( true );
    dihedralAction->setCheckable ( true );
    rotateBondAction->setCheckable(true);
    m_measureActions->setExclusive ( true );


    navigationtoolbar->addAction ( noneAction );
    navigationtoolbar->addAction ( distanceAction );
    navigationtoolbar->addAction ( angleAction );
    navigationtoolbar->addAction ( dihedralAction );
    navigationtoolbar->addAction ( rotateBondAction );
    noneAction->setChecked ( true );

    m_navigationtoolbar.push_front(navigationtoolbar);

    //Add the navigation tool bar to the main window
    for(int i=0;i<m_navigationtoolbar.size();++i)
    {
        addToolBar(Qt::RightToolBarArea, m_navigationtoolbar.at(i));
        m_navigationtoolbar.at(i)->show();
    }
    update();
}


//---------Functions for Open, Read and Update Files:

void KryoMolMainWindow::OnLoadButton()
{
    QString filename = QFileDialog::getOpenFileName ( this );

    if ( filename.isEmpty() )  return;

    OnOpenFile(filename);
}

void KryoMolMainWindow::OnLoadUVFolderAction()
{
    QString foldername=QFileDialog::getExistingDirectory();
    if ( foldername.isEmpty() ) return;
    OpenUVFolder(foldername);

}

void KryoMolMainWindow::OnLoadIRFolderAction()
{
    QString foldername=QFileDialog::getExistingDirectory();
    if ( foldername.isEmpty() ) return;
    OpenIRFolder(foldername);

}

void KryoMolMainWindow::OnOpenFile(QString filename)
{
    m_file=filename;
    OpenFile();
}

void KryoMolMainWindow::OpenFile()
{

    bool b = isGaussianFile();

    if (b)
    {
        QApplication::setOverrideCursor(QCursor(Qt::BusyCursor));
        statusBar()->showMessage(tr("Opening file..."));

        InitGaussian();

        OpenGaussianFile();

        FinishGaussian();

        QApplication::restoreOverrideCursor();
        statusBar()->clearMessage();

    }
    else
    {
        OpenGenericFile();
    }

    this->show();
}


void KryoMolMainWindow::OpenUVFolder(QString foldername)
{
    QDir dir(foldername);
    QFileInfoList flist=dir.entryInfoList(QDir::Files);
    //This should be a list of files containing the ecd computations
    QStringList validfiles;
    bool hasdensity=true;
    bool hasorbitals=true;
    m_world = new kryomol::World(this,World::opticalvisor);
    m_hasdensity=m_hasorbitals=m_hasalphabeta=true;
    for(int i=0;i<flist.size();++i)
    {
        QString f=flist.at(i).absoluteFilePath();
        qDebug() << "file=" << f << endl;
#ifdef __MINGW32__
        kryomol::ParserFactory factory(std::filesystem::u8path(f.toUtf8().data()));
#else
        kryomol::ParserFactory factory(f.toStdString().c_str());
#endif
        kryomol::Parser* qparser=factory.BuildParser();
        if ( !factory.existDensity() )
        {
            m_hasdensity=false;
        }
        if ( !factory.existOrbitals() )
        {
            m_hasorbitals=false;
        }

        if ( !factory.existAlphaBeta() )
        {
            m_hasalphabeta=false;
        }

        for(const auto& j : qparser->Jobs() )
        {
            if ( j.type == kryomol::uv )
            {
                std::vector<Molecule> mol;
                qparser->SetMolecules(&mol);
                qparser->Parse(j.pos);
                qparser->ParseUV(j.pos);
                if ( m_world->Molecules().empty() )
                {
                    m_world->Molecules().push_back(mol.back());
                }
                else
                {
                    m_world->Molecules().back().Frames().push_back(mol.back().Frames().back());
                }

            }
        }
        delete qparser;
    }

    qDebug() << "nframes=" << m_world->Molecules().back().Frames().size() << endl;

    //Should this work ?
    SetBondOrders();
    m_world->Initialize();
    this->InitWidgets(kryomol::uv,m_hasdensity,m_hasorbitals);

}

void KryoMolMainWindow::OpenIRFolder(QString foldername)
{
    QDir dir(foldername);
    QFileInfoList flist=dir.entryInfoList(QDir::Files);
    //This should be a list of files containing the ecd computations
    QStringList validfiles;
    m_world = new kryomol::World(this,World::freqvisor);
    m_hasdensity=m_hasorbitals=m_hasalphabeta=true;
    for(int i=0;i<flist.size();++i)
    {
        QString f=flist.at(i).absoluteFilePath();
#ifdef __MINGW32__
        kryomol::ParserFactory factory(std::filesystem::u8path(f.toUtf8().data()));
#else
        kryomol::ParserFactory factory(f.toStdString().c_str());
#endif
        kryomol::Parser* qparser=factory.BuildParser();
        if ( !factory.existDensity() )
        {
            m_hasdensity=false;
        }
        if ( !factory.existOrbitals() )
        {
            m_hasorbitals=false;
        }

        if ( !factory.existAlphaBeta() )
        {
            m_hasalphabeta=false;
        }

        for(const auto& j : qparser->Jobs() )
        {
            if ( j.type == kryomol::freq )
            {
                std::vector<Molecule> mol;
                qparser->SetMolecules(&mol);
                qparser->Parse(j.pos);
                qparser->ParseFrequencies(j.pos);
                if ( m_world->Molecules().empty() )
                {
                    m_world->Molecules().push_back(mol.back());
                }
                else
                {
                    m_world->Molecules().back().Frames().push_back(mol.back().Frames().back());
                }

            }
        }
        delete qparser;
    }

    qDebug() << "nframes=" << m_world->Molecules().back().Frames().size() << endl;

    //Should this work ?
    SetBondOrders();
    m_world->Initialize();
    this->InitWidgets(kryomol::freq,m_hasdensity,m_hasorbitals);
}

void KryoMolMainWindow::InitWidgets(kryomol::JobType t,bool hasdensity,bool hasorbitals)
{
    //Build the stacked widget for holding the different type of jobs
    delete m_stackedwidget;
    m_stackedwidget = new QStackedWidget (this);
    //If a orca widget was created delete it
    delete m_orcawidget;
    m_orcawidget=nullptr;
    //Get the initial size
    m_size = this->size();


    QJobWidget* q=nullptr;
    if ( t == kryomol::uv )
    {
        q= new QJobUVWidget(m_file, m_world, this);
    }
    else if ( t == kryomol::freq )
    {
        q= new QJobFreqWidget(m_file, m_world, this);
    }

    m_measures = new QMeasureWidget(this);
    connect(m_world, SIGNAL ( currentFrame ( size_t ) ), this, SLOT ( OnUpdateMeasures ( size_t )) );
    connect(m_world->Visor(),SIGNAL(distance(QString&)),m_measures,SLOT(OnWriteDistance(QString&)));
    connect(m_world->Visor(),SIGNAL ( angle ( QString& ) ), m_measures, SLOT ( OnWriteAngle ( QString& ) ) );
    connect(m_world->Visor(),SIGNAL ( dihedral ( QString& ) ), m_measures, SLOT ( OnWriteDihedral ( QString& ) ) );
    connect(m_world, SIGNAL ( currentFrame ( size_t ) ), this, SLOT ( OnUpdateMeasures ( size_t )) );

    connect(m_measures,SIGNAL(clearAll()),m_world->Visor(), SLOT(OnClearMeasures()));
    connect(m_measures,SIGNAL(distanceChange(int)),m_world->Visor(),SLOT(OnDistanceChange(int)));
    connect(m_measures,SIGNAL(angleChange(int)),m_world->Visor(),SLOT(OnAngleChange(int)));
    connect(m_measures,SIGNAL(dihedralChange(int)),m_world->Visor(),SLOT(OnDihedralChange(int)));
    connect(m_measures, SIGNAL(showDistances(bool)),m_world->Visor(),SLOT(OnShowDistances(bool)));
    q->addWidget(m_measures);
    m_measures->hide();
    if ( dynamic_cast<QJobUVWidget*>(q) )
    {
    QUVWidget* uvwidget = ( static_cast<QUVWidget*> (q->widget(1)));
    if ( hasdensity )
    {
        m_orbitals = new QOrbitalWidget(this);
        m_orbitals->SetRenderOrbitals(RenderOrbitals(m_world->Molecules().back().Frames().back()));
        connect(m_world, SIGNAL ( currentFrame ( size_t ) ), this, SLOT ( OnUpdateFrameForOrbitals ( size_t )) );



        connect(m_orbitals,SIGNAL(drawDensity(bool)),this,SLOT(OnDrawDensity(bool)));
        connect(this,SIGNAL(showDensity(bool)),m_world->Visor(),SLOT(OnShowDensity(bool)));

        if ( hasorbitals )
        {
            m_orbitals->SetBeta(m_hasalphabeta);
            uvwidget->SetBeta(m_hasalphabeta);
            m_orbitals->ListOrbitals();

            uvwidget->SetCheckableTransitionChanges();

            connect(m_orbitals,SIGNAL(offtransitions(bool)),uvwidget,SLOT(OffShowTransitionChanges(bool)));
            connect(m_orbitals,SIGNAL(transparenceChange(float)),m_world->Visor(),SLOT(OnTransparenceChange(float)));
            connect(uvwidget,SIGNAL(showtransition(int)),m_orbitals,SLOT(OnShowTransitionChange(int)));
            connect(uvwidget,SIGNAL(showdensities(int)),m_orbitals,SLOT(OnShowDensityChange(int)));
            connect(uvwidget,SIGNAL(offshowtransitions(bool)),m_orbitals,SLOT(OffButtons(bool)));
        }

        q->addWidget(m_orbitals);

        m_measures->hide();
        m_orbitals->hide();
    }
    else
    {
        if (!m_world->Molecules().back().Frames().back().TransitionChanges().empty())
            uvwidget->SetCheckableTransitionCoefficients();
    }
    }

    QList<int> sizelist = q->sizes();

    for(int i=0;i<sizelist.size();++i)
    {
        if ( i == 0 )
        {
        sizelist[0]=m_size.width()-m_size.width()/4;
        }
        else if ( i == 1 )
        {
        sizelist[1]=m_size.width()/4;
        }
        else sizelist[i]=0;
    }
    q->setSizes(sizelist);

    m_stackedwidget->addWidget(q);
    OnLastFrame();
    //Show the first job and update m_world and m_measures
    QJobWidget* jq = ( static_cast<QJobWidget*> ( m_stackedwidget->currentWidget()));
    m_world = jq->GetWorld();
    m_measures = ( static_cast<QMeasureWidget*> ( jq->findChild<QMeasureWidget*>() ));
    m_orbitals = ( static_cast<QOrbitalWidget*> ( jq->findChild<QOrbitalWidget*>() ));
    setCentralWidget(m_stackedwidget);
    FinishGaussian();
    this->show();

}

void KryoMolMainWindow::OnUpdateFrameForOrbitals(size_t f)
{
    m_orbitals->SetRenderOrbitals(RenderOrbitals(m_world->Molecules().back().Frames().at(f)));
    this->update();
}

void KryoMolMainWindow::OnOpenRecentFile()
{
    QAction *action = qobject_cast<QAction *> ( sender() );
    if ( action )
    {
        m_file=action->data().toString();

        bool b = isGaussianFile();

        if (b)
        {
            QApplication::setOverrideCursor(QCursor(Qt::BusyCursor));
            statusBar()->showMessage(tr("Opening file..."));

            InitGaussian();

            OpenGaussianFile();

            FinishGaussian();

            QApplication::restoreOverrideCursor();
            statusBar()->clearMessage();

        }
        else
        {

            OpenGenericFile();
        }

        this->show();
    }
}


void KryoMolMainWindow::OpenGenericFile()
{
    kryomol::KryoMolApplication::SetFile ( m_file );

#ifdef __MINGW32__
    kryomol::ParserFactory factory(std::filesystem::u8path(m_file.toUtf8().data()));
#else
    kryomol::ParserFactory factory(m_file.toStdString().c_str());
#endif


    m_qparser=factory.BuildParser();

    if ( m_qparser == nullptr )
    {
        QMessageBox k;
        k.setText ( "Could not determine file type" );
        k.exec();
        return;
    }

    try
    {
        Parse();

    }
    catch ( ... )
    {
        //first of all clear the world, otherwise the rendering egine will crash on Apply Transformation
        m_world->Clear();
        QMessageBox::warning( this,"","Error during parsing" );
        return;
    }

    SetBondOrders();
    m_world->Initialize();


    if (m_hasdensity)
    {
        m_orbitals->SetRenderOrbitals(RenderOrbitals(m_world->Molecules().back().Frames().back()));

        connect(m_orbitals,SIGNAL(drawDensity(bool)),this,SLOT(OnDrawDensity(bool)));
        connect(m_orbitals,SIGNAL(transparenceChange(float)),m_world->Visor(),SLOT(OnTransparenceChange(float)));
        connect(this,SIGNAL(showDensity(bool)),m_world->Visor(),SLOT(OnShowDensity(bool)));

        if (m_hasorbitals)
        {
            m_orbitals->SetBeta(m_hasalphabeta);
            m_orbitals->ListOrbitals();
        }
    }

    QString apptitle = "KryoMol";
    setWindowTitle ( apptitle + " file://"+m_file );
    UpdateRecentFiles();

    m_mdcontrol->SetWorld ( m_world );
    m_pdbcontrol->SetWorld( m_world );
    m_mdcontrol->Init();
    m_pdbcontrol->Init();

    //Why we have to force update ?
    this->show();

}


void KryoMolMainWindow::OpenGaussianFile()
{

    kryomol::KryoMolApplication::SetFile ( m_file );
#ifdef __MINGW32__
    kryomol::ParserFactory factory(std::filesystem::u8path(m_file.toUtf8().data()));
#else
    kryomol::ParserFactory factory(m_file.toStdString().c_str());
#endif

    m_qparser=factory.BuildParser();

    if ( m_qparser == nullptr )
    {
        QMessageBox k;
        k.setText ( "Could not determine file type" );
        k.exec();
        return;
    }

    m_hasdensity=factory.existDensity();
    m_hasorbitals=factory.existOrbitals();
    m_hasalphabeta=factory.existAlphaBeta();

    //Looking for an Optimization Job:

    QString apptitle = "KryoMol";
    setWindowTitle ( apptitle + " file://"+m_file );

    std::vector<kryomol::JobHeader> jobs=m_qparser->Jobs();

    if ( jobs.empty() ) return;

    for ( std::vector<kryomol::JobHeader>::iterator it=jobs.begin();it!=jobs.end();++it )
    {
        if (it->type == kryomol::singlepoint)
        {
            m_world = new kryomol::World (this,kryomol::World::kryovisor);
            m_qparser->SetMolecules(&m_world->Molecules() );
            try
            {
                m_qparser->Parse(it->pos);
            }
            catch(...)
            {
                std::cerr << "error opening gaussian file";
                delete m_world;
                return;
            }
            m_world->Initialize();
            SetBondOrders();
            QJobSpWidget* q = new QJobSpWidget(m_world,this);
            m_measures = new QMeasureWidget(this);
            m_orbitals = new QOrbitalWidget(this);

            connect(m_world->Visor(),SIGNAL(distance(QString&)),m_measures,SLOT(OnWriteDistance(QString&)));
            connect(m_world->Visor(), SIGNAL ( angle ( QString& ) ), m_measures, SLOT ( OnWriteAngle ( QString& ) ) );
            connect(m_world->Visor(), SIGNAL ( dihedral ( QString& ) ), m_measures, SLOT ( OnWriteDihedral ( QString& ) ) );
            connect(m_world, SIGNAL ( currentFrame ( size_t ) ), this, SLOT ( OnUpdateMeasures ( size_t )) );
            connect(m_measures,SIGNAL(clearAll()),m_world->Visor(), SLOT(OnClearMeasures()));
            connect(m_measures,SIGNAL(distanceChange(int)),m_world->Visor(),SLOT(OnDistanceChange(int)));
            connect(m_measures,SIGNAL(angleChange(int)),m_world->Visor(),SLOT(OnAngleChange(int)));
            connect(m_measures,SIGNAL(dihedralChange(int)),m_world->Visor(),SLOT(OnDihedralChange(int)));
            connect(m_measures, SIGNAL(showDistances(bool)),m_world->Visor(),SLOT(OnShowDistances(bool)));

            if ( m_world->Molecules().back().Frames().size() == 1)
            {
                if ( m_hasdensity )
                {
                    m_orbitals->SetRenderOrbitals(RenderOrbitals(m_world->Molecules().back().CurrentFrame()));

                    connect(m_orbitals,SIGNAL(drawDensity(bool)),this,SLOT(OnDrawDensity(bool)));
                    connect(m_orbitals,SIGNAL(transparenceChange(float)),m_world->Visor(),SLOT(OnTransparenceChange(float)));
                    connect(this,SIGNAL(showDensity(bool)),m_world->Visor(),SLOT(OnShowDensity(bool)));

                    if ( m_hasorbitals )
                    {
                        m_orbitals->SetBeta(m_hasalphabeta);
                        m_orbitals->ListOrbitals();
                    }

                }
            }
            else
            {
                if ( m_hasdensity && !m_hasorbitals )
                {
                    if (!m_world->Molecules().back().CurrentFrame().ElectronicDensityData().Density().Empty())
                    {
                        m_orbitals->SetRenderOrbitals(RenderOrbitals(m_world->Molecules().back().CurrentFrame()));

                        connect(m_orbitals,SIGNAL(drawDensity(bool)),this,SLOT(OnDrawDensity(bool)));
                        connect(m_orbitals,SIGNAL(transparenceChange(float)),m_world->Visor(),SLOT(OnTransparenceChange(float)));
                        connect(this,SIGNAL(showDensity(bool)),m_world->Visor(),SLOT(OnShowDensity(bool)));
                    }
                    else
                    {
                        m_orbitals->HideAllButtons();
                    }

                    connect(m_world, SIGNAL ( currentFrame ( size_t ) ), this, SLOT ( OnUpdateDensities ( size_t )) );

                }

                if ( m_hasorbitals )
                {
                    if (!m_world->Molecules().back().CurrentFrame().OrbitalsData().BasisCenters().empty())
                    {
                        m_orbitals->SetRenderOrbitals(RenderOrbitals(m_world->Molecules().back().CurrentFrame()));

                        connect(m_orbitals,SIGNAL(drawDensity(bool)),this,SLOT(OnDrawDensity(bool)));
                        connect(m_orbitals,SIGNAL(transparenceChange(float)),m_world->Visor(),SLOT(OnTransparenceChange(float)));
                        connect(this,SIGNAL(showDensity(bool)),m_world->Visor(),SLOT(OnShowDensity(bool)));

                        m_orbitals->SetBeta(m_hasalphabeta);
                        m_orbitals->ListOrbitals();
                    }
                    else
                    {
                        m_orbitals->HideAllButtons();
                    }
                    connect(m_world, SIGNAL ( currentFrame ( size_t ) ), this, SLOT ( OnUpdateOrbitals ( size_t )) );
                }

            }


            q->addWidget(m_measures);
            q->addWidget(m_orbitals);
            QList<int> list = q->sizes();
            list[0]=m_size.width()-m_size.width()/4;
            list[1]=m_size.width()/4;

            q->setSizes(list);
            m_measures->hide();
            m_orbitals->hide();
            m_stackedwidget->addWidget(q);
        }

        if (it->type == kryomol::opt || it->type == kryomol::dyn )
        {
            m_world = new World (this,World::optvisor);
            m_qparser->SetMolecules(&m_world->Molecules());
            try
            {
                m_qparser->Parse(it->pos);
            }
            catch(...)
            {
                std::cerr << "error opening gaussian file";
                delete m_world;
                return;
            }
            m_world->Initialize();
            SetBondOrders();
            QJobWidget* q;
            if ( it->type == kryomol::opt )
                q = new QJobOptWidget(m_world,this);
            else q = new QJobDynWidget(m_world,this);

            m_measures = new QMeasureWidget(this);
            m_orbitals = new QOrbitalWidget(this);

            connect(m_world->Visor(),SIGNAL(distance(QString&)),m_measures,SLOT(OnWriteDistance(QString&)));
            connect(m_world->Visor(), SIGNAL ( angle ( QString& ) ), m_measures, SLOT ( OnWriteAngle ( QString& ) ) );
            connect(m_world->Visor(), SIGNAL ( dihedral ( QString& ) ), m_measures, SLOT ( OnWriteDihedral ( QString& ) ) );
            connect(m_world, SIGNAL ( currentFrame ( size_t ) ), this, SLOT ( OnUpdateMeasures ( size_t )) );
            connect(m_measures,SIGNAL(clearAll()),m_world->Visor(), SLOT(OnClearMeasures()));
            connect(m_measures,SIGNAL(distanceChange(int)),m_world->Visor(),SLOT(OnDistanceChange(int)));
            connect(m_measures,SIGNAL(angleChange(int)),m_world->Visor(),SLOT(OnAngleChange(int)));
            connect(m_measures,SIGNAL(dihedralChange(int)),m_world->Visor(),SLOT(OnDihedralChange(int)));
            connect(m_measures, SIGNAL(showDistances(bool)),m_world->Visor(),SLOT(OnShowDistances(bool)));

            if ( m_world->Molecules().back().Frames().size() == 1)
            {
                if ( m_hasdensity )
                {
                    m_orbitals->SetRenderOrbitals(RenderOrbitals(m_world->Molecules().back().Frames().back()));

                    connect(m_orbitals,SIGNAL(drawDensity(bool)),this,SLOT(OnDrawDensity(bool)));
                    connect(m_orbitals,SIGNAL(transparenceChange(float)),m_world->Visor(),SLOT(OnTransparenceChange(float)));
                    connect(this,SIGNAL(showDensity(bool)),m_world->Visor(),SLOT(OnShowDensity(bool)));

                    if ( m_hasorbitals )
                    {
                        m_orbitals->SetBeta(m_hasalphabeta);
                        m_orbitals->ListOrbitals();
                    }

                }
            }
            else
            {
                if ( m_hasdensity && !m_hasorbitals )
                {
                    if (!m_world->Molecules().back().CurrentFrame().ElectronicDensityData().Density().Empty())
                    {
                        m_orbitals->SetRenderOrbitals(RenderOrbitals(m_world->Molecules().back().CurrentFrame()));

                        connect(m_orbitals,SIGNAL(drawDensity(bool)),this,SLOT(OnDrawDensity(bool)));
                        connect(m_orbitals,SIGNAL(transparenceChange(float)),m_world->Visor(),SLOT(OnTransparenceChange(float)));
                        connect(this,SIGNAL(showDensity(bool)),m_world->Visor(),SLOT(OnShowDensity(bool)));
                    }
                    else
                    {
                        m_orbitals->HideAllButtons();
                    }

                    connect(m_world, SIGNAL ( currentFrame ( size_t ) ), this, SLOT ( OnUpdateDensities ( size_t )) );

                }

                if (m_hasorbitals)
                {
                    if (!m_world->Molecules().back().CurrentFrame().OrbitalsData().BasisCenters().empty())
                    {
                        m_orbitals->SetRenderOrbitals(RenderOrbitals(m_world->Molecules().back().CurrentFrame()));

                        connect(m_orbitals,SIGNAL(drawDensity(bool)),this,SLOT(OnDrawDensity(bool)));
                        connect(m_orbitals,SIGNAL(transparenceChange(float)),m_world->Visor(),SLOT(OnTransparenceChange(float)));
                        connect(this,SIGNAL(showDensity(bool)),m_world->Visor(),SLOT(OnShowDensity(bool)));

                        m_orbitals->SetBeta(m_hasalphabeta);
                        m_orbitals->ListOrbitals();
                    }
                    else
                    {
                        m_orbitals->HideAllButtons();
                    }
                    connect(m_world, SIGNAL ( currentFrame ( size_t ) ), this, SLOT ( OnUpdateOrbitals ( size_t )) );
                }

            }


            /*if (existDensity())
               {
                   m_orbitals->SetRenderOrbitals(RenderOrbitals(m_world->Molecules().back().Frames().back()));

                   connect(m_orbitals,SIGNAL(drawDensity(bool)),this,SLOT(OnDrawDensity(bool)));
                   connect(m_orbitals,SIGNAL(transparenceChange(float)),m_world->Visor(),SLOT(OnTransparenceChange(float)));
                   connect(this,SIGNAL(showDensity(bool)),m_world->Visor(),SLOT(OnShowDensity(bool)));

                   if (existOrbitals())
                   {
                       m_orbitals->SetBeta(existAlphaBeta());
                       m_orbitals->ListOrbitals();
                   }
               }*/

            q->addWidget(m_measures);
            m_measures->hide();
            q->addWidget(m_orbitals);
            m_orbitals->hide();
            QList<int> list = q->sizes();
            if ( it->type == kryomol::opt )
            {
                list[0]=m_size.width()-m_size.width()/4;
                list[1]=m_size.width()/4;
                list[2]=0;
                list[3]=0;
            }
            else { //dyn
                list[0]=m_size.width()-m_size.width()/4;
                list[1]=m_size.width()/4;
            }
            q->setSizes(list);
            m_stackedwidget->addWidget(q);
            OnLastFrame();

        }

        if (it->type == kryomol::freq)
        {
            m_world = new kryomol::World(this,World::freqvisor);
            m_qparser->SetMolecules(&m_world->Molecules());
            try
            {
                m_qparser->Parse(it->pos);
            }
            catch( ...)
            {
                std::cerr << "Error parsing gaussian file" << std::endl;
                delete m_world;
                return;
            }

            try
            {
                m_qparser->ParseFrequencies ( it->pos );
            }
            catch(...)
            {
                std::cerr << "Error parsing frequencies" << std::endl;
                return;

            }
            m_world->Initialize();
            SetBondOrders();

            QJobFreqWidget* q = new QJobFreqWidget (m_file, m_world, this);
            m_measures = new QMeasureWidget(this);
            m_orbitals = new QOrbitalWidget(this);

            connect(m_world->Visor(),SIGNAL(distance(QString&)),m_measures,SLOT(OnWriteDistance(QString&)));
            connect(m_world->Visor(),SIGNAL(angle(QString&)),m_measures,SLOT(OnWriteAngle(QString&)));
            connect(m_world->Visor(),SIGNAL(dihedral(QString&)),m_measures,SLOT(OnWriteDihedral(QString&)));
            connect(m_world,SIGNAL(currentFrame(size_t)),this,SLOT(OnUpdateMeasures(size_t)));

            connect(m_measures,SIGNAL(clearAll()),m_world->Visor(),SLOT(OnClearMeasures()));
            connect(m_measures,SIGNAL(distanceChange(int)),m_world->Visor(),SLOT(OnDistanceChange(int)));
            connect(m_measures,SIGNAL(angleChange(int)),m_world->Visor(),SLOT(OnAngleChange(int)));
            connect(m_measures,SIGNAL(dihedralChange(int)),m_world->Visor(),SLOT(OnDihedralChange(int)));
            connect(m_measures, SIGNAL(showDistances(bool)),m_world->Visor(),SLOT(OnShowDistances(bool)));

            if ( m_world->Molecules().back().Frames().size() == 1)
            {
                if ( m_hasdensity )
                {
                    m_orbitals->SetRenderOrbitals(RenderOrbitals(m_world->Molecules().back().Frames().back()));

                    connect(m_orbitals,SIGNAL(drawDensity(bool)),this,SLOT(OnDrawDensity(bool)));
                    connect(m_orbitals,SIGNAL(transparenceChange(float)),m_world->Visor(),SLOT(OnTransparenceChange(float)));
                    connect(this,SIGNAL(showDensity(bool)),m_world->Visor(),SLOT(OnShowDensity(bool)));

                    if ( m_hasorbitals )
                    {
                        m_orbitals->SetBeta(m_hasalphabeta);
                        m_orbitals->ListOrbitals();
                    }

                }
            }
            else
            {
                if ( m_hasdensity && !m_hasorbitals )
                {
                    if (!m_world->Molecules().back().Frames().back().ElectronicDensityData().Density().Empty())
                    {
                        m_orbitals->SetRenderOrbitals(RenderOrbitals(m_world->Molecules().back().CurrentFrame()));

                        connect(m_orbitals,SIGNAL(drawDensity(bool)),this,SLOT(OnDrawDensity(bool)));
                        connect(m_orbitals,SIGNAL(transparenceChange(float)),m_world->Visor(),SLOT(OnTransparenceChange(float)));
                        connect(this,SIGNAL(showDensity(bool)),m_world->Visor(),SLOT(OnShowDensity(bool)));
                    }
                    else
                    {
                        m_orbitals->HideAllButtons();
                    }

                    connect(m_world, SIGNAL ( currentFrame ( size_t ) ), this, SLOT ( OnUpdateDensities ( size_t )) );

                }

                if ( m_hasorbitals)
                {
                    if (!m_world->Molecules().back().Frames().back().OrbitalsData().BasisCenters().empty())
                    {
                        m_orbitals->SetRenderOrbitals(RenderOrbitals(m_world->Molecules().back().CurrentFrame()));

                        connect(m_orbitals,SIGNAL(drawDensity(bool)),this,SLOT(OnDrawDensity(bool)));
                        connect(m_orbitals,SIGNAL(transparenceChange(float)),m_world->Visor(),SLOT(OnTransparenceChange(float)));
                        connect(this,SIGNAL(showDensity(bool)),m_world->Visor(),SLOT(OnShowDensity(bool)));

                        m_orbitals->SetBeta(m_hasalphabeta);
                        m_orbitals->ListOrbitals();
                    }
                    else
                    {
                        m_orbitals->HideAllButtons();
                    }
                    connect(m_world, SIGNAL ( currentFrame ( size_t ) ), this, SLOT ( OnUpdateOrbitals ( size_t )) );
                }

            }

            if ( m_hasdensity )
            {
                m_orbitals->SetRenderOrbitals(RenderOrbitals(m_world->Molecules().back().Frames().back()));

                connect(m_orbitals,SIGNAL(drawDensity(bool)),this,SLOT(OnDrawDensity(bool)));
                connect(m_orbitals,SIGNAL(transparenceChange(float)),m_world->Visor(),SLOT(OnTransparenceChange(float)));
                connect(this,SIGNAL(showDensity(bool)),m_world->Visor(),SLOT(OnShowDensity(bool)));

                if (m_hasorbitals)
                {
                    m_orbitals->SetBeta(m_hasalphabeta);
                    m_orbitals->ListOrbitals();
                }
            }


            q->addWidget(m_measures);
            q->addWidget(m_orbitals);
            QList<int> list = q->sizes();
            list[0]=m_size.width()-m_size.width()/4;
            list[1]=m_size.width()/4;
            list[2]=0;
            list[3]=0;
            q->setSizes(list);
            m_measures->hide();
            m_orbitals->hide();
            m_stackedwidget->addWidget(q);
            OnLastFrame();
        }

        if (it->type == kryomol::uv)
        {
            m_world = new kryomol::World(this,World::opticalvisor);
            m_qparser->SetMolecules(&m_world->Molecules());
            try
            {
                m_qparser->Parse(it->pos);
            }
            catch( ...)
            {
                std::cerr << "Error parsing gaussian file" << std::endl;
                delete m_world;
                return;
            }

            m_world->Initialize();
            SetBondOrders();
            m_qparser->ParseUV ( it->pos );

            QJobUVWidget* q = new QJobUVWidget(m_file, m_world, this);

            m_measures = new QMeasureWidget(this);
            m_orbitals = new QOrbitalWidget(this);

            connect(m_world->Visor(),SIGNAL(distance(QString&)),m_measures,SLOT(OnWriteDistance(QString&)));
            connect(m_world->Visor(),SIGNAL ( angle ( QString& ) ), m_measures, SLOT ( OnWriteAngle ( QString& ) ) );
            connect(m_world->Visor(),SIGNAL ( dihedral ( QString& ) ), m_measures, SLOT ( OnWriteDihedral ( QString& ) ) );
            connect(m_world, SIGNAL ( currentFrame ( size_t ) ), this, SLOT ( OnUpdateMeasures ( size_t )) );

            connect(m_measures,SIGNAL(clearAll()),m_world->Visor(), SLOT(OnClearMeasures()));
            connect(m_measures,SIGNAL(distanceChange(int)),m_world->Visor(),SLOT(OnDistanceChange(int)));
            connect(m_measures,SIGNAL(angleChange(int)),m_world->Visor(),SLOT(OnAngleChange(int)));
            connect(m_measures,SIGNAL(dihedralChange(int)),m_world->Visor(),SLOT(OnDihedralChange(int)));
            connect(m_measures, SIGNAL(showDistances(bool)),m_world->Visor(),SLOT(OnShowDistances(bool)));

            QUVWidget* uvwidget = ( static_cast<QUVWidget*> (q->widget(1)));
            if ( m_hasdensity )
            {
                m_orbitals->SetRenderOrbitals(RenderOrbitals(m_world->Molecules().back().Frames().back()));

                connect(m_orbitals,SIGNAL(drawDensity(bool)),this,SLOT(OnDrawDensity(bool)));
                connect(this,SIGNAL(showDensity(bool)),m_world->Visor(),SLOT(OnShowDensity(bool)));

                if ( m_hasorbitals )
                {
                    m_orbitals->SetBeta(m_hasalphabeta);
                    uvwidget->SetBeta(m_hasalphabeta);
                    m_orbitals->ListOrbitals();

                    uvwidget->SetCheckableTransitionChanges();

                    connect(m_orbitals,SIGNAL(offtransitions(bool)),uvwidget,SLOT(OffShowTransitionChanges(bool)));
                    connect(m_orbitals,SIGNAL(transparenceChange(float)),m_world->Visor(),SLOT(OnTransparenceChange(float)));
                    connect(uvwidget,SIGNAL(showtransition(int)),m_orbitals,SLOT(OnShowTransitionChange(int)));
                    connect(uvwidget,SIGNAL(showdensities(int)),m_orbitals,SLOT(OnShowDensityChange(int)));
                    connect(uvwidget,SIGNAL(offshowtransitions(bool)),m_orbitals,SLOT(OffButtons(bool)));
                }
            }
            else
            {
                if (!m_world->Molecules().back().Frames().back().TransitionChanges().empty())
                    uvwidget->SetCheckableTransitionCoefficients();
            }

            q->addWidget(m_measures);
            q->addWidget(m_orbitals);
            QList<int> list = q->sizes();
            list[0]=m_size.width()-m_size.width()/4;
            list[1]=m_size.width()/4;
            list[2]=0;
            list[3]=0;
            q->setSizes(list);
            m_measures->hide();
            m_orbitals->hide();
            m_stackedwidget->addWidget(q);
            OnLastFrame();
        }
    }


    //Show the first job and update m_world and m_measures
    QJobWidget* q = ( static_cast<QJobWidget*> ( m_stackedwidget->currentWidget()));
    m_world = q->GetWorld();
    m_measures = ( static_cast<QMeasureWidget*> ( q->findChild<QMeasureWidget*>() ));
    m_orbitals = ( static_cast<QOrbitalWidget*> ( q->findChild<QOrbitalWidget*>() ));
    setCentralWidget(m_stackedwidget);

}


void KryoMolMainWindow::Parse()
{
    m_world->Clear();
    m_qparser->SetMolecules ( &m_world->Molecules() );
    m_qparser->Parse();
}

bool KryoMolMainWindow::isGaussianFile()
{
    kryomol::KryoMolApplication::SetFile ( m_file );

#ifdef __MINGW32__
    kryomol::ParserFactory factory(std::filesystem::u8path(m_file.toUtf8().data()));
#else
    kryomol::ParserFactory factory(m_file.toStdString().c_str());
#endif

    bool ftype = factory.isGaussianFile();

    return ftype;
}

bool KryoMolMainWindow::isGaussianFile(QString file)
{

#ifdef __MINGW32__
    kryomol::ParserFactory factory(std::filesystem::u8path(m_file.toUtf8().data()));
#else
    kryomol::ParserFactory factory(m_file.toStdString().c_str());
#endif

    return factory.isGaussianFile();
}

void KryoMolMainWindow::UpdateRecentFiles()
{

    QSettings settings;
    QStringList files = settings.value ( "RecentFiles" ).toStringList();
    if ( !m_file.isEmpty() )
    {
        files.removeAll ( m_file );
        files.prepend ( m_file );
        while ( files.size() > maxrecentfiles )
            files.removeLast();
        settings.setValue ( "RecentFiles", files );
    }
    int nfiles = qMin ( files.size(), maxrecentfiles );

    for ( int i = 0; i < nfiles; ++i )
    {
        QString text = tr ( "&%1 %2" ).arg ( i + 1 ).arg ( QFileInfo ( files[i] ).fileName() );
        m_recentfileactions.at ( i )->setText ( text );
        m_recentfileactions.at ( i )->setData ( files[i] );
        m_recentfileactions.at ( i )->setVisible ( true );
    }
    for ( int j = nfiles; j < maxrecentfiles; ++j )
        m_recentfileactions[j]->setVisible ( false );

}

void KryoMolMainWindow::SetBondOrders()
{
    for ( std::vector<kryomol::Molecule>::iterator it=m_world->Molecules().begin();it!=m_world->Molecules().end();++it )
    {
        it->SetBonds();
    }
}



//-------------- Functions for Selecting Frame:

void KryoMolMainWindow::OnFirstFrame()
{
    m_world->SelectFrame ( 0 );
}

void KryoMolMainWindow::OnLastFrame()
{
    if ( m_world->CurrentMolecule() )
        m_world->SelectFrame ( m_world->CurrentMolecule()->Frames().size()-1 );
}

void KryoMolMainWindow::OnNextFrame()
{
    if ( m_world->CurrentMolecule() )
    {
        if ( m_world->CurrentMolecule()->CurrentFrameIndex() < m_world->CurrentMolecule()->Frames().size() )
        {
            m_world->SelectFrame ( m_world->CurrentMolecule()->CurrentFrameIndex() +1 );
        }
    }
}

void KryoMolMainWindow::OnPreviousFrame()
{
    //take care of unsigned size_t
    if ( m_world->CurrentMolecule() )
    {
        if ( m_world->CurrentMolecule()->CurrentFrameIndex() > 0 )
        {
            m_world->SelectFrame ( m_world->CurrentMolecule()->CurrentFrameIndex()-1 );
        }
    }
}

//-------------- Functions for Measuring Distance, Angles and Dihedrals:

void KryoMolMainWindow::OnNoneSelection()
{
    m_world->Visor()->OnNoneSelection();

}

void KryoMolMainWindow::OnMeasureDistances()
{
    m_world->Visor()->OnMeasureDistances();

}

void KryoMolMainWindow::OnMeasureAngles()
{
    m_world->Visor()->OnMeasureAngles();

}

void KryoMolMainWindow::OnMeasureDihedrals()
{
    m_world->Visor()->OnMeasureDihedrals();

}

void KryoMolMainWindow::OnRotateBond()
{
    m_world->Visor()->OnRotateBond();

}

void KryoMolMainWindow::OnShowMeasures()
{
    if (m_stackedwidget->currentWidget()->findChild<QJobOptWidget*>() ||  m_stackedwidget->currentWidget()->findChild<QJobFreqWidget*>())
    {

        QJobWidget* q = ( static_cast<QJobWidget*> ( m_stackedwidget->currentWidget()));
        QList<int> list = q->sizes();
        if ( m_measures->isHidden() )
        {
            list[list.size()-2] = m_size.width()/6;
            list[0] = list[0] - m_size.width()/6;
            q->setSizes(list);
            m_measures->show();
        }
        else
        {
            list[0] = list[0] + m_size.width()/6;
            list[list.size()-2] = 0;
            q->setSizes(list);
            m_measures->hide();
        }
    }
    else
    {
        QSplitter* q = ( static_cast<QSplitter*>(m_stackedwidget->currentWidget()));
        QList<int> list = q->sizes();
        if ( m_measures->isHidden() )
        {
            list[list.size()-2] = m_size.width()/6;
            list[0] = list[0] - m_size.width()/6;
            q->setSizes(list);
            m_measures->show();
        }
        else
        {
            list[0] = list[0] + m_size.width()/6;
            list[list.size()-2] = 0;
            q->setSizes(list);
            m_measures->hide();
        }
    }
}

void KryoMolMainWindow::OnUpdateMeasures(size_t frame)
{
    QStringList listDistances = m_world->Visor()->GetDistances(frame);
    QStringList listAngles = m_world->Visor()->GetAngles(frame);
    QStringList listDihedrals = m_world->Visor()->GetDihedrals(frame);

    if (!listDistances.empty())
        m_measures->updateDistances(listDistances);
    if (!listAngles.empty())
        m_measures->updateAngles(listAngles);
    if (!listDihedrals.empty())
        m_measures->updateDihedrals(listDihedrals);

}

////-------------- Functions for Represent Orbital Information:


void KryoMolMainWindow::OnDrawDensity(bool b)
{
    if (b)
    {
        m_world->CurrentMolecule()->CurrentFrame().SetPositiveDensity(m_orbitals->GetDensity().PositiveRenderDensityVector());
        m_world->CurrentMolecule()->CurrentFrame().SetNegativeDensity(m_orbitals->GetDensity().NegativeRenderDensityVector());

    }
    emit showDensity(b);
}

void KryoMolMainWindow::OnUpdateDensities(size_t frame)
{
    if (!m_world->Molecules().back().CurrentFrame().ElectronicDensityData().Density().Empty())
        m_orbitals->SetRenderOrbitals(RenderOrbitals(m_world->Molecules().back().CurrentFrame()));
    else
        m_orbitals->HideAllButtons();


}

void KryoMolMainWindow::OnUpdateOrbitals(size_t frame)
{
    if (!m_world->Molecules().back().CurrentFrame().OrbitalsData().BasisCenters().empty())
    {
        m_orbitals->SetRenderOrbitals(RenderOrbitals(m_world->Molecules().back().CurrentFrame()));
        m_orbitals->SetBeta(m_hasalphabeta);
        m_orbitals->ListOrbitals();
    }
    else
        m_orbitals->HideAllButtons();

}


// ------------- Functions for Export Information:

void KryoMolMainWindow::OnExportVectorGraphics()
{
    QString filter1= "Encasulated Postscript (*.eps)";
    QString filter2= "Postscript (*.ps)";
    QString filter3= "SVG (*.svg)";
    QString filter4="PDF (*.pdf)";
    QStringList filters;
    QString filter;
    filters << filter1 << filter2 << filter3 << filter4;
    QString s=QFileDialog::getSaveFileName ( this,"Save Vector Graphics", QDir::home().canonicalPath(), filter1+";;"+filter2+";;"+filter3+";;"+filter4,&filter );


    if ( s.isEmpty() )
        return;


    kryomol::StringTokenizer path ( s.toStdString(),"." );
    std::string back=kryomol::toupper ( path.back() );
    QString mode=QString ( path.back().c_str() );
    if ( ( back !="EPS" ) && ( back !="PS" ) && ( back !="SVG" ) )
    {
        if ( filter==filter1 )
        {
            mode="EPS";
            s+=".eps";
        }
        if ( filter==filter2 )
        {
            mode="PS";
            s+=".ps";
        }
        if ( filter==filter3 )
        {
            mode="SVG";
            s+=".svg";
        }

    }

    QApplication::setOverrideCursor ( QCursor ( Qt::WaitCursor ) );
    m_world->Visor()->ExportVectorPicture ( s,mode );
    QApplication::restoreOverrideCursor();
}

void KryoMolMainWindow::OnExportRasterGraphics()
{

    QList<QByteArray> formats= QImageWriter::supportedImageFormats();
    QString filters;
    QList<QByteArray>::iterator it;
    for ( it=formats.begin();it!=formats.end()-1;++it )
    {
        filters+=QString ( *it );
        filters+=";;";
    }
    filters+=formats.back();

    QString filter;
    QString s=QFileDialog::getSaveFileName ( this,"Save Raster Graphics", QDir::home().canonicalPath(), filters ,&filter );

    if ( s.isEmpty() ) return;

    kryomol::StringTokenizer path ( s.toStdString(),"." );
    std::string b=kryomol::toupper ( path.back() );
    std::string fst=filter.toStdString();
    if ( b != fst )
    {
        s = s +"." +filter;
    }


    QImage img=m_world->Visor()->grabFrameBuffer();
    QFile* f= new QFile ( s );
    QImageWriter w ( f,filter.toStdString().c_str() );
    if ( ! w.write ( img ) )
    {
        QMessageBox b;
        b.setText ( w.errorString() );
        b.exec();
    }
    delete f;
}


void KryoMolMainWindow::OnExportGeomCurrent()
{
    if ( m_world->CurrentMolecule() )
    {
        std::ostringstream   f;
        kryomol::XYZWriter w(m_world->CurrentMolecule(),false);
        f << w;
        QClipboard* cp=QApplication::clipboard();
        if ( cp )
        {
            cp->setText(QString(f.str().c_str()));
        }

    }
}

void KryoMolMainWindow::OnExportGeom()
{
    if ( m_world->CurrentMolecule() )
    {
        std::ostringstream   f;
        kryomol::XYZWriter w(m_world->CurrentMolecule(),true);
        f << w;
        QClipboard* cp=QApplication::clipboard();
        if ( cp )
        {
            cp->setText(QString(f.str().c_str()));
        }

    }

}

void KryoMolMainWindow::OnExportACES()
{

    QClipboard* c=QApplication::clipboard();

    AcesWriter awr ( m_world->CurrentMolecule() );
    std::stringstream s;

    s << awr;
    c->setText ( s.str().c_str() );

}

void KryoMolMainWindow::OnExportCPMD()
{
    QClipboard* c=QApplication::clipboard();
    CPMDWriter w (m_world->CurrentMolecule() );
    std::stringstream s;
    s << w;
    c->setText ( s.str().c_str() );
}


//void MobDevWindow::OnExportGaussianNICSGrid()
//{
/*   std::pair<Coordinate,Coordinate> box= m_world->CurrentMolecule()->CurrentFrame().Box();


   QNICSGridDialog d(box,this);
   if ( d.exec() == QDialog::OK )
   {
       box=d.Box();

       QClipboard* c=QApplication::clipboard();
       if ( c == NULL )
       {
         QMessageBox::foregroundRole();
         return;
       }
       std::stringstream s;
       s << "#Route geom=connectivity" << "\n" << "title" << std::endl;
       s << "0 1" << std::endl;
       size_t natoms=m_world->CurrentMolecule()->Atoms().size();
       const std::vector<Atom>& atoms=m_world->CurrentMolecule()->Atoms();
       const std::vector<Coordinate& c>=m_world->CurrentMolecule()->CurrentFrame().XYZ();
       //write the geometry
       for ( size_t i=0:i<natoms;++i )
       {

         s << std::setiosflags ( std::ios::fixed ) << std::setw ( 3 ) << setiosflags ( std::ios::left ) << it->Symbol()
         << std::setw ( 10 )  << std::resetiosflags ( std::ios::left )
         << std::setiosflags ( std::ios::right ) << std::setiosflags ( std::ios::fixed ) << std::setprecision ( 6 )
         << ct->x() <<
         std::setw ( 10 ) <<ct->y()
         << std::setw ( 10 ) << ct->z()
         << std::resetiosflags ( std::ios::right ) << std::endl;
       }
       //generate now the connectivity
       int width=static_cast<int>(floor(std::log10(natoms)));
       for(size_t i=0;i<natoms;++i)
       {
           s << std::setiosflags ( std::ios::fixed ) << setiosflags ( std::ios::left );
           s << std::setw(width) << (i+1) << std::endl;
           std::vector<size_t> bonds=m_world->CurrentMolecule()->Bonds().at(i);
           for(std::vector::const_iterator bt=bonds.begin();bt!=bonds.end();++bt)
           {
               s << std::setw(width) << (*bt+1) << std::setw(3) << 1;
           }
           s << std::endl;
       }

   }
*/

//}
void KryoMolMainWindow::OnExportGaussian ( bool withhessian )
{
    QClipboard* c=QApplication::clipboard();
    if ( c == NULL )
    {
        std::cerr << "Could not open clipboard";
    }

    std::stringstream s;
    s << 0 <<"," << 1 << std::endl;
    if ( !withhessian )
    {
        std::cout << "Export Gaussian geometry" << std::endl;

        s << m_world->CurrentMolecule()->CurrentFrame();
        s << std::endl;

        c->setText ( s.str().c_str() );

        return;
    }

    std::cout << "Export Gaussian geometry and forces" << std::endl;

    std::vector<Atom>::iterator it;
    std::vector<Coordinate>::const_iterator ct=m_world->CurrentMolecule()->InputOrientation().begin();
    for ( it=m_world->CurrentMolecule()->Atoms().begin();it!=m_world->CurrentMolecule()->Atoms().end();++it,++ct )
    {

        s << std::setiosflags ( std::ios::fixed ) << std::setw ( 3 ) << std::setiosflags ( std::ios::left ) << it->Symbol()
          << std::setw ( 10 )  << std::resetiosflags ( std::ios::left )
          << std::setiosflags ( std::ios::right ) << std::setiosflags ( std::ios::fixed ) << std::setprecision ( 6 )
          << ct->x() <<
             std::setw ( 10 ) <<ct->y()
          << std::setw ( 10 ) << ct->z()
          << std::resetiosflags ( std::ios::right ) << std::endl;
    }
    s << std::endl << std::endl;
    s << std::resetiosflags ( std::ios::fixed );

    s << std::setiosflags ( std::ios::fixed );
    s << std::setiosflags ( std::ios::uppercase ) << std::setprecision ( 8 );


    D1Array<double>& forces=m_world->CurrentMolecule()->CurrentFrame().GetForces();
    D2Array<double>& hessian=m_world->CurrentMolecule()->CurrentFrame().GetHessian();

    if ( forces.size() )
    {
        size_t i=0;
        int j;
        while ( i < forces.size() )
        {
            for ( j=0;j<6;++j )
            {
                s << std::setw ( 12 )  << forces ( i++ );
                if ( i == forces.size() )
                    break;
            }
            s << std::endl;

        }

    }

    if ( hessian.NRows() )
    {
        size_t i,j;
        short k=0;
        for ( i=0;i<hessian.NRows();i++ )
            for ( j=0;j<=i;j++ )
            {
                s << std::setw ( 12 ) << hessian ( i,j );
                if ( ++k > 5 )
                {
                    s << std::endl;
                    k=0;
                }
            }

    }
    s << std::endl;
    c->setText ( s.str().c_str() );

}


void KryoMolMainWindow::OnExportGaussianForces(bool withhessian)
{
    OnExportGaussian(true);
}


void KryoMolMainWindow::OnSetGraphMode(QAction* act)
{

    m_world->Visor()->SetGraphMode(static_cast<kryomol::GLVisorBase::graphmode>(act->data().toInt()));
}

void KryoMolMainWindow::OnGraphModeChanged(kryomol::GLVisorBase::graphmode mode)
{
    for(int waction=0;waction<m_gmodegroup->actions().size();++waction)
    {
        if ( m_gmodegroup->actions().at(waction)->data().toInt() == static_cast<int>(mode) )
        {
            m_gmodegroup->actions().at(waction)->setChecked(true);
            return;
        }
    }
}

void KryoMolMainWindow::OnCenterVisiblePart()
{
    m_world->Visor()->CenterMolecule(false);
}

void KryoMolMainWindow::OnCenterWholeMolecule()
{
    m_world->Visor()->CenterMolecule(true);

}

void KryoMolMainWindow::OnChangedAnalysisMenu(QAction* act)
{
    int i=act->data().toInt();
    OnChangePlugin(i);
}

void KryoMolMainWindow::OnChangePlugin ( int i )
{
    m_uistack->setCurrentIndex(i);

    //we need to restore state of the dock otherwise layout gets corrupted i dont know why
    //#ifdef Q_WS_MAC
    //  restoreDockWidget(m_plugindock); //this is necessary even in  4.4.2
    //#endif //we need to restore state of the dock otherwise layout gets corrupted i dont know why

    m_currentplugin=i;
    m_world->SetCurrentPlugin ( m_app->Plugins() [m_currentplugin] );
    m_world->Visor()->update();
}

void KryoMolMainWindow::OnAbout()
{
    QString str;
    QTextStream stream(&str);
    stream << m_app->applicationName() << endl;
    stream << "Revision:" << kryomol::KryoMolInfo::Version() << endl;
    stream << "svn:" << kryomol::KryoMolInfo::SvnVersion() << endl;
    stream << endl;
    stream << "Main Author:" << endl;

    QStringList authors=kryomol::KryoMolInfo::Authors();
    for(QList<QString>::iterator it=authors.begin();it!=authors.end();++it)
    {
        stream << (*it) << endl;
    }
    stream << endl;
    stream << "Collaborators:" << endl;

    QStringList coll=kryomol::KryoMolInfo::Collaborators();
    for(QList<QString>::iterator it=coll.begin();it!=coll.end();++it)
    {
        stream << (*it) << endl;
    }
    stream << endl;
    stream << "Kryomol makes use of the following 3d party open-source libraries:";
    stream << endl;
    QStringList tparty=kryomol::KryoMolInfo::ThirdPartyPackages();
    for(QList<QString>::iterator it=tparty.begin();it!=tparty.end();++it)
    {
        stream << (*it) << endl;
    }
    stream << endl;
    stream << kryomol::KryoMolInfo::Copyright() << endl;
    QMessageBox::about(this,"",str);
    stream << endl;
    stream << kryomol::KryoMolInfo::Disclaimer() << endl;
}

//-------------- Functions of Widgets for Generic Files:

void KryoMolMainWindow::OnOpenConfButton()
{
    QSplitter* q = ( static_cast<QSplitter*>(m_stackedwidget->currentWidget()));
    QList<int> list = q->sizes();

    if ( m_mdcontrol->isHidden() )
    {
        list[2] = m_size.width()/4;
        list[0] = list[0] - m_size.width()/4;
        q->setSizes(list);
        m_mdcontrol->show();
    }
    else
    {
        list[0] = list[0] + m_size.width()/4;
        list[2] = 0;
        q->setSizes(list);
        m_mdcontrol->hide();
    }
}


void KryoMolMainWindow::OnOpenPDBButton()
{
    QSplitter* q = ( static_cast<QSplitter*>(m_stackedwidget->currentWidget()));
    QList<int> list = q->sizes();

    if ( m_pdbcontrol->isHidden() )
    {
        list[1] = m_size.width()/6;
        list[0] = list[0] - m_size.width()/6;
        q->setSizes(list);
        m_pdbcontrol->show();
    }
    else
    {
        list[0] = list[0] + m_size.width()/6;
        list[1] = 0;
        q->setSizes(list);
        m_pdbcontrol->hide();
    }
}


void KryoMolMainWindow::OnOpenPluginButton()
{
    if ( m_uistack->isHidden() ) m_uistack->show();
    else m_uistack->hide();
}

void KryoMolMainWindow::OnOpenOrbitalButton(bool b)
{
    QSplitter* q = ( static_cast<QSplitter*>(m_stackedwidget->currentWidget()));
    QList<int> list = q->sizes();

    if ( b )
    {
        list[list.count()-1] = m_size.width()/3;
        list[0] = list[0] - m_size.width()/3;
        q->setSizes(list);
        m_orbitals->show();
    }
    else
    {
        list[0] = list[0] + m_size.width()/3;
        list[list.count()-1] = 4;
        q->setSizes(list);
        m_orbitals->hide();
    }
}


//------------- Functions of Actions for Gaussian Files:

void KryoMolMainWindow::OnStartAnimation()
{
    QJobWidget* current = ( static_cast<QJobWidget*> ( m_stackedwidget->currentWidget()));
    if ( ! current ) return;
    if (KryoVisorOpt* w = current->findChild<KryoVisorOpt*>())
        w->OnStartAnimation();
    if (KryoVisorFreq* w = current->findChild<KryoVisorFreq*>())
        w->OnStartAnimation();
}


void KryoMolMainWindow::OnStopAnimation()
{
    QJobWidget* current = ( static_cast<QJobWidget*> ( m_stackedwidget->currentWidget()));
    if ( ! current ) return;
    if (KryoVisorOpt* w = current->findChild<KryoVisorOpt*>())
        w->OnStopAnimation();
    if (KryoVisorFreq* w = current->findChild<KryoVisorFreq*>())
        w->OnStopAnimation();

}

void KryoMolMainWindow::OnChangeJob()
{
    if ((m_stackedwidget->currentIndex()>-1)&&(m_stackedwidget->count()>1))
    {
        QJobWidget* q;
        if (m_stackedwidget->currentIndex() == m_stackedwidget->count()-1)
        {
            m_stackedwidget->setCurrentIndex(0);
            q = ( static_cast<QJobWidget*> ( m_stackedwidget->currentWidget()));
            m_world = q->GetWorld();
            m_measures = ( static_cast<QMeasureWidget*> ( q->findChild<QMeasureWidget*>() ));
            m_orbitals = ( static_cast<QOrbitalWidget*> ( q->findChild<QOrbitalWidget*>() ));
        }
        else
        {
            m_stackedwidget->setCurrentIndex(m_stackedwidget->currentIndex()+1);
            q = ( static_cast<QJobWidget*> ( m_stackedwidget->currentWidget()));
            m_world = q->GetWorld();
            m_measures = ( static_cast<QMeasureWidget*> ( q->findChild<QMeasureWidget*>() ));
            m_orbitals = ( static_cast<QOrbitalWidget*> ( q->findChild<QOrbitalWidget*>() ));
        }
        setCentralWidget(m_stackedwidget);
    }


    m_measureActions->actions().at(0)->setChecked(true);

    FinishGaussian();
}

void KryoMolMainWindow::OnChangeJobN()
{
    int job;
    for (int i=0; i<m_jobs->actions().count(); i++)
    {
        if (m_jobs->actions().at(i)->isChecked())
            job = i;
    }
    if ((m_stackedwidget->currentIndex()>-1)&&(m_stackedwidget->count()>1))
    {
        QJobWidget* q;

        m_stackedwidget->setCurrentIndex(job);
        q = ( static_cast<QJobWidget*> ( m_stackedwidget->currentWidget()));
        m_world = q->GetWorld();
        m_measures = ( static_cast<QMeasureWidget*> ( q->findChild<QMeasureWidget*>() ));
        m_orbitals = ( static_cast<QOrbitalWidget*> ( q->findChild<QOrbitalWidget*>() ));

        setCentralWidget(m_stackedwidget);
    }


    m_measureActions->actions().at(0)->setChecked(true);

    FinishGaussian();
}

void KryoMolMainWindow::OnProtonateTrigonalCenterAction()
{
    QMessageBox box;
    box.setText("Protonate trigonal center\n"
                "click on central atom and then the three surrounding atoms\
                to add the proton in a clock-wise or counterclockwise way");
                if ( box.exec() == QMessageBox::Ok )
                {
                    connect(m_world->Visor(),SIGNAL(selectedatoms(std::vector<size_t> )),
                    this,SLOT(OnProtonateTrigonalCenter(std::vector<size_t> )));
                    m_world->Visor()->OnFourAtomsSelectionMode();

                }

}

void KryoMolMainWindow::OnProtonateTrigonalCenter(std::vector<size_t> selatoms)
{
    std::vector<Atom>& atoms=m_world->CurrentMolecule()->Atoms();
    std::vector<Coordinate>& c=m_world->CurrentMolecule()->CurrentFrame().XYZ();
    std::vector<Bond>* b=nullptr;
    if ( m_world->CurrentMolecule()->Bonds().empty() )
    {
        b=&m_world->CurrentMolecule()->CurrentFrame().Bonds();
    } else b=&m_world->CurrentMolecule()->Bonds();
    //get the first vector
    qDebug() << "i=" << selatoms[0] << "," << selatoms[1] << "," << selatoms[2];
    Coordinate c1=c[selatoms[1]]-c[selatoms[0]];
    Coordinate c2=c[selatoms[2]]-c[selatoms[0]];
    Coordinate ch=c1^c2;
    double norm=ch.Norm();
    ch/=norm;
    ch+=c[selatoms[0]];
    atoms.push_back(Atom(1));
    c.push_back(ch);
    m_world->Visor()->update();
    Bond hc(selatoms[0],atoms.size()-1);
    b->push_back(hc);
}

void KryoMolMainWindow::OnRunOrcaWidget()
{
    qDebug() << "OnRunOrcaWidget";
    //qDebug() << "orca widget" << m_orcawidget;
    delete m_orcawidget;

    m_orcawidget = new OrcaDialog(this);
    m_orcawidget->SetMolecule(this->m_world->CurrentMolecule());
    connect(m_orcawidget,SIGNAL(outputfile(QString )),this,SLOT(OnOpenFile(QString )));
    m_orcawidget->exec();
}



