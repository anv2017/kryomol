/*****************************************************************************************
                            kryovisor.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <sstream>
#include <iomanip>
#include <algorithm>

#include <qclipboard.h>
#include <qapplication.h>
#include <qmenu.h>
#include <qtimer.h>
//Added by qt3to4:
#include <QMouseEvent>
#include <QWheelEvent>
#include <QInputDialog>
#include <QMessageBox>

#ifdef Q_OS_MAC
#include <OpenGl/glu.h>
#else
#include <GL/glu.h>
#endif

#include "kryovisor.h"
#include "qtimeropenmenu.h"
#include "aceswriter.h"
#include "animation.h"
#include "gl2ps.h"
#include "xyzwriter.h"
#include "qrenumberoptionsdialog.h"
#include "world.h"

using namespace kryomol;

const float rotationpass=5.0f;

class QMenu;

class KryoVisor::KryoVisorPrivate
{
  public:
    KryoVisorPrivate() {}
    ~KryoVisorPrivate() {}
    //draw wireframe when moving
    bool m_bwfonmoving;
};


//---------------------- KryoVisor Functions ---------------------------------------

KryoVisor::KryoVisor ( World* world, QWidget* parent, const QGLWidget* shareWidget, Qt::WindowFlags f )
    :GLVisor (world, parent,shareWidget,f ) , m_world ( world )
{
    _d = new KryoVisorPrivate();
    _d->m_bwfonmoving=true;

    m_selectionmode=NONE;
    m_bshowespcharges=false;

    m_framespeed=100;
    m_screencontrols.resize ( 6 );
    setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Expanding);

}


KryoVisor::~KryoVisor()
{
    delete _d;
}

void KryoVisor::Initialize()
{
  GLVisor::Initialize();
  m_timer= new QTimer ( this );


}


void KryoVisor::SetSelectionMode ( selectionmode mode )
{
  m_selectionmode=mode;
}

void KryoVisor::OnRecalculateConnectivity()
{
  m_world->CurrentMolecule()->SetBonds();
  update();
}

void KryoVisor::OnMoveToCentroid()
{
  m_world->CurrentMolecule()->CalculateMassCenter();
  update();
}
void KryoVisor::OnSelection ( QMouseEvent* e )
{
    GLVisor::OnSelection(e);
}

void KryoVisor::OnResetSelection ( bool repaint /*=true*/ )
{
  m_selatoms.clear();
  m_world->CurrentMolecule()->ResetSelection();
  if ( repaint )
    updateGL();
}


void KryoVisor::RenderScreenText()
{
  std::stringstream label;
  label << "#" << m_world->CurrentMolecule()->CurrentFrameIndex()+1 << " of " << m_world->CurrentMolecule()->Frames().size();
  QColor col ( QColor ( 255,255,255 ) );
  qglColor ( col );
  renderText ( rect().left() +30,rect().bottom() - 30,QString ( label.str().c_str() ),GLFont() );
}

void KryoVisor::OnFirstFrame()
{
  m_world->SelectFrame ( 0 );
  RefreshDistances();
}

void KryoVisor::OnLastFrame()
{
  if ( m_world->CurrentMolecule() )
    m_world->SelectFrame ( m_world->CurrentMolecule()->Frames().size()-1 );
  RefreshDistances();
}

void KryoVisor::OnNextFrame()
{
  if ( m_world->CurrentMolecule() )
  m_world->SelectFrame ( m_world->CurrentMolecule()->CurrentFrameIndex() +1 );
  RefreshDistances();
}

void KryoVisor::OnPreviousFrame()
{
  //take care of unsigned size_t
  if ( m_world->CurrentMolecule() )
  if ( m_world->CurrentMolecule()->CurrentFrameIndex() > 0 )
  {
    m_world->SelectFrame ( m_world->CurrentMolecule()->CurrentFrameIndex()-1 );
  }
  RefreshDistances();
}

void KryoVisor::RefreshDistances()
{
  std::vector<Molecule::pair>::iterator it;

  for ( it=m_rpairs.begin();it!=m_rpairs.end();it++ )
  {
    QString str;
    str.sprintf ( "%.3f",Coordinate::Distance ( m_world->CurrentMolecule()->CurrentFrame().XYZ().at(it->i),m_world->CurrentMolecule()->CurrentFrame().XYZ().at(it->j) ) );
    it->value=str.toStdString();
  }

  update();
}

void KryoVisor::mousePressEvent ( QMouseEvent* e )
{
  GLVisorBase::mousePressEvent ( e );
  switch ( e->button() )
  {
    case Qt::LeftButton:
      OnSelection ( e );
      break;
    case Qt::RightButton:
      break;
    default:
      break;
  }
}


void KryoVisor::wheelEvent ( QWheelEvent* e )
{
    GLVisor::wheelEvent(e);
}


void KryoVisor::mouseMoveEvent ( QMouseEvent* e )
{
    GLVisor::mouseMoveEvent ( e );
}


void KryoVisor::OnMeasureDistances ( bool b )
{
  if ( b )
    SetSelectionMode ( DISTANCES );
  else
    SetSelectionMode ( NONE );
  GLVisor::OnMeasureDistances();
  OnResetSelection();
}


void KryoVisor::Gromacs3dout ( size_t i, size_t j, size_t k, size_t l )
{
  Coordinate c1=m_world->CurrentMolecule()->CurrentFrame().XYZ().at(i);
  Coordinate c2=m_world->CurrentMolecule()->CurrentFrame().XYZ().at(j);
  Coordinate c3=m_world->CurrentMolecule()->CurrentFrame().XYZ().at(k);
  Coordinate c4=m_world->CurrentMolecule()->CurrentFrame().XYZ().at(l);

  Coordinate rij=c2-c1;
  Coordinate rjk=c3-c2;
  Coordinate ril=c4-c1;
  Coordinate rop=c2^c3;

  //output in nanometers
  float a=0.1*Coordinate::ScalarProduct(ril,rij)/rij.Norm();
  float b=0.1*Coordinate::ScalarProduct(ril,rjk)/rjk.Norm();
  float c=0.1*Coordinate::ScalarProduct(ril,rop)/rop.Norm();
  std::stringstream ss;
  ss << i+1 << " " << j+1 << " " << k+1 << " " << l+1 << " " << a << " " << b << " " <<  c << std::endl;
  QString s ( ss.str().c_str() );
  QMessageBox information;
  information.setInformativeText(s);
}


void KryoVisor::OnMeasureAngles ( bool b )
{
  if ( b )
    SetSelectionMode ( ANGLES );
  else
    SetSelectionMode ( NONE );
  GLVisor::OnMeasureAngles();
  OnResetSelection();
}


void KryoVisor::OnMeasureDihedrals ( bool b )
{
  if ( b )
    SetSelectionMode ( DIHEDRALS );
  else
    SetSelectionMode ( NONE );
  GLVisor::OnMeasureDihedrals();
  OnResetSelection();
}


void KryoVisor::OnRotateBond (bool b)
{
  if ( b )
    SetSelectionMode ( ROTATEBOND );
  else
    SetSelectionMode ( NONE );
  GLVisor::OnRotateBond();
  OnResetSelection();
}


void KryoVisor::RenderScene(GLenum mode)
{
    GLVisor::RenderScene(mode);
}


void KryoVisor::RenderMolecule ( size_t index, GLenum mode/*=GL_RENDER*/ )
{
    GLVisor::RenderMolecule ( index,mode );
}

void KryoVisor::OnShowESPCharges()
{

  m_bshowespcharges=!m_chargemenu->actions().at(0)->isChecked();;
  m_chargemenu->actions().at(0)->setChecked(m_bshowespcharges );
  update();
}


void KryoVisor::OnClearDistances()
{
  m_rpairs.clear();
  update();
}


void KryoVisor::OnSelectPoint ( size_t current )
{
  if ( current < 0 || current >= m_world->Molecules().size() )
    return;

  m_world->SelectMolecule(current);
}


void KryoVisor::RenumberAtoms()
{
  if ( m_selectionmode == eRenumberHeavy )
  {
    std::vector<size_t>  hatoms;

    {
      for ( size_t it=0;it!=m_selatoms.size();++it )
      {
        for (size_t bt=0; bt<m_world->CurrentMolecule()->Bonds().size(); bt++)
        {
            if ((m_world->CurrentMolecule()->Bonds().at(bt).I() == it)&(m_world->CurrentMolecule()->Atoms().at(m_world->CurrentMolecule()->Bonds().at(bt).J()).Symbol()=="H"))
                hatoms.push_back(m_world->CurrentMolecule()->Bonds().at(bt).J());
        }
      }
    }
    std::vector<size_t>::const_iterator at;
    for ( at=hatoms.begin();at!=hatoms.end();++at )
    {
      m_selatoms.push_back ( *at );
    }
  }
  size_t counter=0;
  std::vector<Atom> renumberedatoms ( m_selatoms.size() );
  for (size_t it=0;it<m_selatoms.size();++it,counter++ )
  {
    renumberedatoms.at(counter)=m_world->CurrentMolecule()->Atoms().at(it);
  }
  m_world->CurrentMolecule()->Atoms()=renumberedatoms;
  m_world->CurrentMolecule()->SetBonds();
}

void KryoVisor::Center()
{
    GLVisor::Center();
}



//---------------------- KryoVisorOpt Functions ---------------------------------------//


KryoVisorOpt::KryoVisorOpt ( World* world, QWidget* parent, const QGLWidget* shareWidget, Qt::WindowFlags f  ) : KryoVisor ( world, parent, shareWidget, f )
{
  m_forcescale=1;

  m_bshowforces=false;
}


void KryoVisorOpt::Initialize()
{
  KryoVisor::Initialize();

  connect ( m_timer,SIGNAL ( timeout() ),this,SLOT ( OnChangeFrame() ) );

  Molecule& firstmolecule=m_world->Molecules().front();
  float maxforce=0.000001;
  std::vector<Coordinate>::const_iterator it;
  for ( it=firstmolecule.CurrentFrame().Gradient().begin();it!=firstmolecule.CurrentFrame().Gradient().end();++it )
  {
    float force=it->Norm();
    if ( force > maxforce )
      maxforce=force;
  }
  //Maximum cylinder should be about 3 A
  m_forcemagnitude=3/maxforce;
}


void KryoVisorOpt::RenderScreenText()
{
  KryoVisor::RenderScreenText();
}


void KryoVisorOpt::OnShowForces(bool b)
{
  m_bshowforces=b;
  update();
}

void KryoVisorOpt::RenderScene(GLenum mode)
{
    KryoVisor::RenderScene(mode);
}

void KryoVisorOpt::RenderMolecule ( size_t index,GLenum mode )
{
    KryoVisor::RenderMolecule(index,mode);
    Molecule& molecule= m_world->Molecules().at(index);

    std::vector<Coordinate>::const_iterator mit;
    const Frame& frame=molecule.CurrentFrame();


    std::vector<Coordinate>::const_iterator ct=frame.XYZ().begin();
    int i=0;
    float color[4];
    color[0]=1.0;
    color[1]=0.5;
    color[2]=0.0;
    color[3]=0.5;

if ( m_bshowforces ) //dont draw anything for wireframe in render mode
{
  glPushMatrix();
  Handlers() [index].ApplyTransformation(molecule.CurrentFrameIndex());
  Coordinate zaxis ( 0,0,1 );
  GLUquadricObj* quadric= gluNewQuadric();
  for ( mit=frame.Gradient().begin();mit!=frame.Gradient().end();++mit,i++,++ct )
  {
          glPushMatrix();
          Coordinate normal=zaxis^(*mit);
          float angle=Coordinate::Angle ( zaxis,(*mit) );
          glTranslatef ( ct->x(),ct->y(),ct->z() );

          glPushMatrix();
          glRotatef ( angle*180/M_PI,normal.x(),normal.y(),normal.z() );
          const int spheres=40;
          Coordinate cyl=(*mit) *m_forcemagnitude*m_forcescale;
          gluCylinder ( quadric,0.05,0.05,cyl.Norm(),spheres,spheres );
          glPopMatrix();

          glPushMatrix();
          {

            glTranslatef ( cyl.x(),cyl.y(),cyl.z() );


            glRotatef ( angle*180/M_PI,normal.x(),normal.y(),normal.z() );
            glMaterialfv ( GL_FRONT_AND_BACK, GL_AMBIENT, color );
            glMaterialfv ( GL_FRONT_AND_BACK, GL_DIFFUSE, color );
            gluCylinder ( quadric,0.09,0.001,0.2,spheres,spheres );
          }
          glPopMatrix();

          glPopMatrix();
      }
  gluDeleteQuadric(quadric);



 glPopMatrix();

}
return;
}


void KryoVisorOpt::OnChangeFrame()
{
  if ( m_world->CurrentMolecule()->CurrentFrameIndex() == ( m_world->CurrentMolecule()->Frames().size() -1 ) )
  {
    m_timer->stop();
    emit playing ( false );
    return;
  }
  //We are doing deep copy, this should be changed on next version
  OnNextFrame();
#ifdef __GNUC__
#warning performance issue here
#endif

}


void KryoVisorOpt::wheelEvent ( QWheelEvent* e )
{
  KryoVisor::wheelEvent ( e );
}

void KryoVisorOpt::OnStartAnimation()
{
  if ( m_world->CurrentMolecule()->CurrentFrameIndex() == ( m_world->CurrentMolecule()->Frames().size() -1 ) )
  {
    m_world->SelectFrame(0);
    emit selectedPoint ( m_world->CurrentMolecule()->CurrentFrameIndex() );
  }
  m_timer->start ( m_framespeed );
  emit playing ( true );
}


void KryoVisorOpt::OnStopAnimation()
{
  m_timer->stop();
  emit playing ( false );
}

//---------------------- KryoVisorFreq Functions ---------------------------------------


KryoVisorFreq::KryoVisorFreq ( World* world, QWidget* parent, const QGLWidget* shareWidget, Qt::WindowFlags f  ) : KryoVisor ( world, parent, shareWidget, f )
{
  m_screencontrols.resize ( 2 );

}

KryoVisorFreq::~KryoVisorFreq()
{
  m_animation->Stop();
}

void KryoVisorFreq::Initialize()
{
  KryoVisor::Initialize();
  m_animation= new Animation ( *m_world->CurrentMolecule() );
  connect(m_animation,SIGNAL( shot() ),this,SLOT( update() ));
  m_animation->BackupMolecule();

}

void KryoVisorFreq::OnPushImage()
{

}

void KryoVisorFreq::OnDistortion ( int v )
{
  m_animation->SetFrame ( v );
  update();
}

void KryoVisorFreq::OnStartAnimation()
{
  m_animation->Start();
  emit playing ( true );
}

bool KryoVisorFreq::IsPlaying() const
{
    return m_animation->IsActive();
}

void KryoVisorFreq::OnStopAnimation()
{
    m_animation->Stop();
    update();
    emit playing ( false );
}

void KryoVisorFreq::SetMode ( int mode,int frame )
{

  m_animation->SetMode ( mode,frame );
  update();
}

void KryoVisorFreq::RenderScreenText()
{
  std::stringstream label;
  label << "Current active vibrational mode: #" << m_animation->GetMode() +1;
  QString str ( label.str().c_str() );
  QColor col ( QColor ( 255,255,255 ) );
  qglColor ( col );
  renderText ( rect().left() +30,rect().bottom() - 30,str,GLFont() );
}

void KryoVisorFreq::RenderScene(GLenum mode)
{
    KryoVisor::RenderScene(mode);
}

void KryoVisorFreq::RenderMolecule ( size_t index, GLenum mode )
{
    KryoVisor::RenderMolecule ( index, mode );
}

void KryoVisorFreq::wheelEvent ( QWheelEvent* e )
{
  KryoVisor::wheelEvent ( e );
}
