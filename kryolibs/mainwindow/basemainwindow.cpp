/*****************************************************************************************
                            basemainwindow.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <sstream>
#include <stdlib.h>
#include <iostream>


#include <qpixmap.h>
#include <QKeyEvent>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QTabWidget>

#include <qlayout.h>
#include <qobject.h>
#include <qimage.h>
#include <qpushbutton.h>
#include <qclipboard.h>
#include <qlibrary.h>
#include <qmessagebox.h>
#include <qapplication.h>
#include <qcombobox.h>
#include <QAction>
#include <qmenubar.h>
#include <qsplitter.h>
#include <qtooltip.h>

#include <QStackedWidget>
#include <QWhatsThis>
#include <QFileDialog>
#include <QSettings>
#include <QUrl>
#include <QTextStream>
#include <QImageWriter>
#include <QStatusBar>


#include "xyzwriter.h"
#include "parserfactory.h"
#include "basemainwindow.h"
#include "qryomolapp.h"
#include "world.h"
#include "plugin.h"
#include "stringtools.h"
#include "parsers.h"
#include "atom.h"
#include "qmolecularlistcontrol.h"
#include "qcwtcombobox.h"
#include "url.h"
#include "qhiddendockwindow.h"
#include "qryomolapp.h"
#include "glvisor.h"
#include "qpdbcontrol.h"
#include "scripter.h"
#include "qorbitalwidget.h"

using namespace kryomol;

const int maxrecentfiles=5;

BaseMainWindow::BaseMainWindow ( QWidget* parent, Qt::WindowFlags f ) : QMainWindow ( parent,f )
{
    setAttribute(Qt::WA_DeleteOnClose);
    InitControls();
}


BaseMainWindow::BaseMainWindow ( const char* file, QWidget* parent, Qt::WindowFlags f ) : QMainWindow ( parent,f ) , m_file ( file )
{
  InitControls();
}

void BaseMainWindow::InitControls()
{
   m_mdcontrol = new QMolecularListControl(this);
   m_pdbcontrol= new QPDBControl(this);
   m_uistack = new QStackedWidget(this);
}

BaseMainWindow::~BaseMainWindow()
{}


kryomol::Plugin* BaseMainWindow::CurrentPlugin()
{
  if ( m_app->Plugins().empty() )
    return NULL;

  return m_app->Plugins().at ( m_currentplugin );
}


void BaseMainWindow::LoadPlugins()
{
  m_app->LoadPlugins();
  if ( m_app->Plugins().empty() )
  {
    QMessageBox err(this);
    err.setText ( "Any plugin seems to be installed on your system" );
    err.exec();
    return;
  }
  std::vector<Plugin*>::iterator mt;
  int i=0;
  m_pluginactions = new QActionGroup(this);


  for ( mt=m_app->Plugins().begin();mt!=m_app->Plugins().end();++mt,++i )
  {
#ifndef Q_WS_MAC
    m_pluginscombo->addItem ( QString ( ( *mt )->Name().c_str() ) );//QString(pJ->HTMLDescription().c_str()));
#endif
    ( *mt )->SetWorld ( m_world );
    m_uistack->addWidget ( (*mt )->UI());
	
    QAction* action=m_analysismenu->addAction(QString ( ( *mt )->Name().c_str() ) );
    action->setCheckable(true);
    action->setData(QVariant(i));
    m_pluginactions->addAction( action );
    ( *mt )->SetMainWindow ( this );
  }
  
  connect(  m_analysismenu,SIGNAL(triggered(QAction* )),this,SLOT(OnPluginTriggered(QAction* )));

  OnChangedPluginsList ( 0 );

}

void BaseMainWindow::keyPressEvent ( QKeyEvent* e )
{
  switch ( e->key() )
  {
    case Qt::Key_F1:
      if ( e->modifiers() == Qt::ShiftModifier )
        QWhatsThis::enterWhatsThisMode();
      break;
    case Qt::Key_V:
      if ( e->modifiers() == Qt::ControlModifier )
        OnPaste();
      break;
    default:
      break;
  }


  QMainWindow::keyPressEvent ( e );

}


void BaseMainWindow::OpenFile()
{
  KryoMolApplication::SetFile ( m_file );
  kryomol::ParserFactory factory(m_file.toStdString().c_str());

  m_qparser=factory.BuildParser();

    if ( m_qparser == NULL )
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

  //Initialize the plugins
  std::vector<Plugin*>::iterator it;
  for ( it=m_app->Plugins().begin();it!=m_app->Plugins().end();++it )
  {
    std::cerr << ( *it )->Name() << std::endl;
    ( *it )->Initialize();
  }
  setWindowTitle ( m_app->applicationName() + " file://"+m_file );
  UpdateRecentFiles();

  m_mdcontrol->SetWorld ( m_world );
  m_pdbcontrol->SetWorld( m_world );
  m_mdcontrol->Init();
  m_pdbcontrol->Init();

  /*if ( m_world->CurrentMolecule()->Frames().size() > 1 )
  {
    m_mldock->show();
  }*/


  /*if ( ! m_world->CurrentMolecule()->Residues().empty() )
  {
    m_pdbdock->show();
  }*/


}

//void BaseMainWindow::OpenFileGaussian()
//{

// NMRDevApplication::SetFile ( m_file );
// qryomol::ParserFactory factory(m_file.toStdString().c_str());

// m_qparser=factory.BuildParser();

//   if ( m_qparser == NULL )
//   {
//     QMessageBox k;
//     k.setText ( "Could not determine file type" );
//     k.exec();
//     return;
//   }


// try
// {
//   ParseGaussian();
//   std::vector<KryoVisor*>::iterator it;
//   for ( it=m_visors.begin();it!=m_visors.end();++it )
//   {
//     ( *it )->SetSelectionMode ( bkmode );
//   }

// }

// catch ( ... )
// {
// //first of all clear the world, otherwise the rendering egine will crash on Apply Transformation
//   m_world->Clear();
//   QMessageBox::warning( this,"","Error during parsing" );
//   return;
// }


//}

//void BaseMainWindow::ParseGaussian()
//{

//  if ( m_qparser )
//  {
//      std::vector<JobHeader> jobs=m_qparser->GetJobs();

//    if ( jobs.empty() ) return; //oops some error

//    std::vector<JobHeader>::iterator it;
//    for ( it=jobs.begin();it!=jobs.end();++it )
//    {
//      //Determine the type of job and proceed
//      switch ( it->type )
//      {
//        case JobType::opt:
//          m_convcontrol = new QConvWidget ( this);
//          break;
//        case JobType::freq:
//#warning change to getURL?
//          m_freqcontrol = new QFreqWidget ( localFilePath(),this );
//          //this is important, otherwise the spectrum gets very small
//          connect ( w,SIGNAL ( mode ( int,int ) ),visor,SLOT ( SetMode ( int,int ) ) );
//          connect ( w,SIGNAL ( movie ( bool ) ),visor,SLOT ( RecordMovie ( bool ) ) );
//          connect ( w,SIGNAL ( distort ( int ) ),visor,SLOT ( OnDistortion ( int ) ) );
//          connect ( w,SIGNAL ( reset() ),visor,SLOT ( OnStopAnimation() ) );
//          connect ( w,SIGNAL ( Type ( IRSpectrum::SpectrumType ) ),this,SLOT ( OnIRTypeChanged ( IRSpectrum::SpectrumType ) ) );

//          break;
//        default:
//          w=NULL;
//      }
//      connect ( visor,SIGNAL ( text ( const QString& ) ),this, SIGNAL ( setStatusBarText ( const QString& ) ) );
//      connect ( visor,SIGNAL ( playing ( bool ) ),this,SLOT ( OnUpdateControls ( bool ) ) );
//      m_visors.push_back ( visor );
//      m_qparser->SetMolecule ( visor->GetMolecule() );

//      visor->setSizePolicy ( QSizePolicy::Expanding,QSizePolicy::Expanding );
//      connect ( visor,SIGNAL ( text ( const QString& ) ),this,SLOT ( OnInfoText ( const QString& ) ) );
//      if ( w ) w->setSizePolicy ( QSizePolicy::Preferred,QSizePolicy::Expanding );

//      if ( m_qparser->ParseFile ( it->pos ) )
//      {
//        switch ( it->type )
//        {
//          case JobType::opt:
//          {
//            QConvWidget* cw=static_cast<QConvWidget*> ( w );
//            cw->SetNData ( visor->GetMolecule()->size() );
//            double* energies=cw->GetEnergies();
//            double* s2= cw->GetS2();
//            double *rmsforce=cw->GetRMSForces();
//            double* maximumforce=cw->GetMaximumForces();
//            double* rmsdisplacement=cw->GetRMSDisplacements();
//            double* maximumdisplacement=cw->GetMaximumDisplacements();
//            cw->SetEnergyLevel ( visor->GetMolecule()->front().EnergyLevel().c_str() );

//            if ( energies )
//            {
//              std::vector<Molecule>::iterator mt;
//              int i=0;
//              mt=visor->GetMolecule()->begin();
//              cw->SetThreshold ( mt->GetThreshold() );
//              for ( i=0;mt!=visor->GetMolecule()->end();mt++,i++ )
//              {
//                energies[i]=mt->GetEnergy();
//                s2[i]=mt->GetS2();
//                rmsforce[i]=mt->GetRMSForce();
//                maximumforce[i]=mt->GetMaximumForce();
//                rmsdisplacement[i]=mt->GetRMSDisplacement();
//                maximumdisplacement[i]=mt->GetMaximumDisplacement();
//              }

//            }
//            visor->Initialize();
//            cw->SetupCurves();
//            connect ( cw,SIGNAL ( selectedPoint ( int ) ),visor,SLOT ( OnSelectPoint ( int ) ) );
//            connect ( visor,SIGNAL ( selectedPoint ( int ) ),cw,SLOT ( OnSelectedPoint ( int ) ) );
//            cw->OnSelectedPoint ( visor->GetCurrentPoint() );
//            connect ( cw,SIGNAL ( forcescale ( float ) ),visor,SLOT ( OnForceScale ( float ) ) );
//            EnableNavigationActions( visor->GetMolecule()->size() > 1 );
//          }
//          break;
//          case JobType::freq:
//          {

//            try
//            {
//              m_qparser->ParseFrequencies ( it->pos );
//              QFreqWidget* cw=static_cast<QFreqWidget*> ( w );
//              connect ( cw,SIGNAL ( showspectrum ( bool ) ),this,SLOT ( OnShowSpectrum ( bool ) ) );

//              visor->Initialize();
//              cw->SetFrequencies ( visor->m_molecule.GetFrequencies() );
//              cw->InitTable();

//            }
//            catch ( ... ) //allow to get opt jobs before broken freq jobs
//            {
//              std::cerr << "Error parsing frequencies" << std::endl;
//            }
//            EnableNavigationActions(true);
//          }
//          break;
//          default:
//            visor->Initialize();

//        }



//      }

//    }
//    //m_kryotab->addTab ( currentTab , m_url.fileName() );
//    //m_kryotab->setTabToolTip ( currentTab, m_url.prettyUrl() );
//    m_tab->setCurrentPage ( 0 );
//    if ( m_bisinsidehtml )
//    {
//      if ( m_tab->count() == 1 )
//        m_tab->hideTabBar();
//    }


//    QList<KryoVisor*> childlist =m_tab->currentPage()->findChildren<KryoVisor*>();

//    if ( !childlist.isEmpty() )
//    {
//      kDebug() << "m_currentvisor=" << m_currentvisor << endl;
//      m_currentvisor=dynamic_cast<KryoVisor*> ( childlist.first() );
//    }


//  }

//  /* The m_bcheckglsettings flag is stored into the konqueror session (KPart is loaded in memory, so it wont be displayed when navigating back and forth*/
//  if ( ! m_bcheckglsettings )
//  {
//    CheckGLSettings();
//    m_bcheckglsettings=true;
//  }

//  return true;
//}



void BaseMainWindow::UpdateRecentFiles()
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

void BaseMainWindow::OnOpenRecentFile()
{
  QAction *action = qobject_cast<QAction *> ( sender() );
  if ( action )
  {
    m_file=action->data().toString();

    OpenFile();
  }

}

bool BaseMainWindow::isGaussianFile()
{
  KryoMolApplication::SetFile ( m_file );
  kryomol::ParserFactory factory(m_file.toStdString().c_str());

  bool ftype = factory.isGaussianFile();

  return ftype;

}

bool BaseMainWindow::existDensity()
{
  KryoMolApplication::SetFile ( m_file );
  kryomol::ParserFactory factory(m_file.toStdString().c_str());

  bool ftype = factory.existDensity();

  return ftype;

}

void BaseMainWindow::Parse()
{
  m_world->Clear();
  m_qparser->SetMolecules ( &m_world->Molecules() );
  m_qparser->Parse();

}

void BaseMainWindow::SetBondOrders()
{
  std::vector<Molecule>::iterator it;
  for ( it=m_world->Molecules().begin();it!=m_world->Molecules().end();++it )
  {
    it->SetBonds();
  }
}


void BaseMainWindow::OnChangePlugin ( int i )
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

void BaseMainWindow::OnPrintButton()
{
//#warning to implement

}


void BaseMainWindow::OnNMRDevOnline()
{
  OpenURL( m_app->URL() );
}

void BaseMainWindow::OnHelp()
{
  QString spldir;

#ifdef Q_WS_WIN

  spldir = qApp->applicationDirPath();
  spldir.truncate ( spldir.lastIndexOf ( "bin" ) );
  spldir +="doc";

#endif

#ifdef Q_OS_UNIX
  char* nppath=getenv ( "KRYOMOL_DIR" );

  QString basedir ( nppath );

  spldir=basedir +"/doc";
  std::cout << spldir.toStdString() << std::endl;
#endif

  spldir +="KryoMol.pdf";

  OpenURL ( spldir );
  
}

void BaseMainWindow::DropHandler ( const QMimeData* src )
{
  qDebug() << "reimplement";
//  QString textfile;
//  //OK lets suppose we have a file list
//  if ( src->hasUrls() )
//  {
//    QList<QUrl> filelist=src->urls();
//    if ( !filelist.isEmpty() )
//    {
//      m_file=filelist.front().toLocalFile();
//      OpenFile();
//      return;
//    }
//  }
//  else
//  {

//    //Try first to get the archive from ChemDraw
//    textfile=QString ( src->data ( "Gaussian Input" ) );
//    if ( textfile.isEmpty() )
//      textfile =QString ( src->text() );

//  }


//  if ( textfile.isEmpty() )
//  {
//    QMessageBox ok;
//    ok.setText ( tr ( "No data in clipboard" ) );
//    ok.exec();
//  }
//  else
//  {
//    //std::stringstream* pf=new std::stringstream(std::ios::binary);
//    std::stringstream f ( std::ios::in|std::ios::out|std::ios::binary );

//    f << textfile.toStdString();
//    std::cout << textfile.toStdString() << "\n" << "Pasted" << "\n" << f.str() << std::endl;



//    qryomol::ParserFactory factory(&f);

//    m_qparser=factory.BuildParser();

//    if ( m_qparser == NULL )
//    {
//      QMessageBox bwar;
//      bwar.setText ( tr ( "No valid data in clipboard" ) );
//      bwar.exec();
//      return;
//    }

//    try
//    {
//      Parse();

//    }
//    catch ( ... )
//    {
//      m_world->Clear();
//      QMessageBox::warning(this,""," Error during parsing" );
//      return;
//    }
//    SetBondOrders();
//    m_world->Initialize();

//    //Initialize the plugins
//    std::vector<Plugin*>::iterator it;
//    for ( it=m_app->Plugins().begin();it!=m_app->Plugins().end();++it )
//    {
//      std::cerr << ( *it )->Name() << std::endl;
//      ( *it )->Initialize();
//    }
//    setWindowTitle ( m_app->applicationName() + " From clipboard" );
//    m_mdcontrol->SetWorld ( m_world );
//    m_mdcontrol->Init();
//    if ( m_world->CurrentMolecule()->Frames().size() > 1 )
//    {
//      m_mdcontrol->show();
//    }


//  }
}
void BaseMainWindow::OnPaste()
{
  const QMimeData* src= QApplication::clipboard()->mimeData();
  if ( src )
    DropHandler ( src );

}

void BaseMainWindow::OnAbout()
{
  QString str;
  QTextStream stream(&str);
  stream << m_app->applicationName() << endl;
  stream << m_app->Version() << endl;
  stream << "Compiled with Qt version" << QT_VERSION_STR << endl;
  stream << "Author " << m_app->Authors() << endl;
  QMessageBox::about(this,"",str);
}


void BaseMainWindow::OnChangedPluginsList ( int i )
{

  for ( int j=0;j<m_pluginactions->actions().size();++j )
  {
    if ( i == j ) m_pluginactions->actions().at(i)->setChecked(true);
  }
  OnChangePlugin ( i );
}

void BaseMainWindow::OnChangedAnalysisMenu(QAction* act)
{
  int i=act->data().toInt();
  OnChangePlugin(i);
}

void BaseMainWindow::OnPluginTriggered(QAction* act)
{
#ifndef Q_WS_MAC
  int i=act->data().toInt();
  m_pluginscombo->setCurrentIndex ( i );
#endif
  
}

void BaseMainWindow::OnExportVectorGraphics()
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


  StringTokenizer path ( s.toStdString(),"." );
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

void BaseMainWindow::OnExportRasterGraphics()
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

  StringTokenizer path ( s.toStdString(),"." );
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

void BaseMainWindow::OnFirstFrame()
{
  m_world->SelectFrame ( 0 );
}

void BaseMainWindow::OnLastFrame()
{
  if ( m_world->CurrentMolecule() )
    m_world->SelectFrame ( m_world->CurrentMolecule()->Frames().size()-1 );
}

void BaseMainWindow::OnNextFrame()
{
  if ( m_world->CurrentMolecule() )
  m_world->SelectFrame ( m_world->CurrentMolecule()->CurrentFrameIndex() +1 );
}

void BaseMainWindow::OnPreviousFrame()
{
  //take care of unsigned size_t
  if ( m_world->CurrentMolecule() )
  if ( m_world->CurrentMolecule()->CurrentFrameIndex() > 0 )
  {
    m_world->SelectFrame ( m_world->CurrentMolecule()->CurrentFrameIndex()-1 );
  }
}

void BaseMainWindow::OnCenterVisiblePart()
{
  m_world->Visor()->CenterMolecule(false);
}

void BaseMainWindow::OnCenterWholeMolecule()
{
  m_world->Visor()->CenterMolecule(true);
  
}

void BaseMainWindow::OnExportGeom()
{
  if ( m_world->CurrentMolecule() )
 {
    std::ostringstream   f;
    XYZWriter w(m_world->CurrentMolecule(),false);
   f << w;
   QClipboard* cp=QApplication::clipboard();
   if ( cp )
  {
    cp->setText(QString(f.str().c_str()));
  }
  
 }
}

void BaseMainWindow::OnSetGraphMode(QAction* act)
{

  m_world->Visor()->SetGraphMode(static_cast<GLVisorBase::graphmode>(act->data().toInt()));
}

void BaseMainWindow::OnGraphModeChanged(kryomol::GLVisorBase::graphmode mode)
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


void BaseMainWindow::OnShowCommandLine(bool b)
{
    if (b)
    {
        if ( m_scripter->CommandLine()->parent() == static_cast<QObject*>(statusBar()) )
            m_scripter->CommandLine()->show();
        else
            statusBar()->addPermanentWidget(m_scripter->CommandLine());
    }
    else m_scripter->CommandLine()->hide();
}



