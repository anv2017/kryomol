/*****************************************************************************************
                            glvisorbase.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <iomanip>

#include "glvisorbase.h"
#include "qapplication.h"
#include <qmessagebox.h>
#include <qcursor.h>
#include <qimage.h>
#include <qclipboard.h>
//Added by qt3to4:
#include <QKeyEvent>
#include <QPixmap>
#include <QMouseEvent>
#include <QWheelEvent>
#include <QColorDialog>
#include <QPainter>
#include <QTimer>

#include <qtimeropenmenu.h>

#include "thermo.h"
#include "gl2ps.h"

#ifdef Q_OS_MAC
#include <OpenGl/glu.h>
#else
#include <GL/glu.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


class kryomol::GLVisorBasePrivate
{
  public:
    GLVisorBasePrivate() :  m_spheres ( 8 )
    {
      m_lgrabbed=false;
      m_rgrabbed=false;
      m_mgrabbed=false;
      m_bismoving=false;
      m_bkcolor=QColor ( 0,0,0 );
      m_glfont.setStyleStrategy ( QFont::OpenGLCompatible );
    }

    ~GLVisorBasePrivate()
    {}

    bool m_vectorgraphics;
    bool m_bshowsymbols;
    bool m_bshownumbers;
    bool m_bshowpdbinfo;
    bool m_bshowdipole;
    bool m_bshowcell;
    bool m_bshowaxis;
    bool m_bshowdensity;
    QFont m_glfont;
    int m_spheres;
    int m_cylinders;
    bool m_rgrabbed;
    bool m_lgrabbed;
    bool m_mgrabbed;
    bool m_bismoving;
    QColor m_bkcolor;
    GLVisorBase::graphmode m_graphmode;
    Coordinate m_camera;
    std::vector<MoleculeHandler> m_handlers;
    Coordinate m_rotationvector;
    std::vector<Molecule>* m_molecules;
    QPoint m_oldpoint;
    GLfloat m_aspect;
    size_t m_currentmolecule;

};

using namespace kryomol;

/**Constructor*/
GLVisorBase::GLVisorBase ( QWidget * parent, const QGLWidget* shareWidget, Qt::WindowFlags f ) :
    QGLWidget ( parent,shareWidget, f ) , m_parent ( parent )
{
  Init();
}

/** Destructor*/
void GLVisorBase::Init()
{
  _d = new GLVisorBasePrivate();
  _d->m_vectorgraphics=false;
  _d->m_currentmolecule=0;
  m_range=5.0f;

  setFocusPolicy ( Qt::ClickFocus );

  _d->m_graphmode=STICKS;
  _d->m_bshowsymbols=false;
  _d->m_bshownumbers=false;
  _d->m_bshowpdbinfo=false;
  _d->m_bshowdipole=false;
  _d->m_bshowcell=false;
  _d->m_bshowaxis=false;
  _d->m_bshowdensity=false;

}

/** destructor*/
GLVisorBase::~GLVisorBase()
{
  delete _d;
}

/** \brief initialize rendering properties

Initialize the lights and blending properties
*/
void GLVisorBase::initializeGL()
{
  glEnable ( GL_NORMALIZE );
  GLfloat light_position[]={ 0.8,0.7,1.0,0.0 };

  GLfloat ambient_color[]={ 0.20,0.20,0.20,1.0 };
  GLfloat specular_color[]= { 1.0,1.0,1.0,1.0 };
  GLfloat difusse_color[]= { 1.0,1.0,1.0,1.0 };

  qglClearColor ( _d->m_bkcolor );

  glLightfv ( GL_LIGHT0, GL_POSITION, light_position );
  glLightfv ( GL_LIGHT0, GL_AMBIENT, ambient_color );
  glLightfv ( GL_LIGHT0, GL_SPECULAR, specular_color );
  glLightfv ( GL_LIGHT0, GL_DIFFUSE, difusse_color );



  glEnable ( GL_LIGHTING );   // so lighting models are used
  glEnable ( GL_LIGHT0 );  // we'll use LIGHT0

  glEnable ( GL_DEPTH_TEST ); // allow z-buffer display

  //blending stuff
  glEnable ( GL_BLEND );   // enable alpha-channel blending
  glBlendFunc ( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

}

/** \brief resize the widget
     Resize the widgedt to width \a w and height \a h    
*/
void GLVisorBase::resizeGL ( int w, int h )
{

  GLsizei width,height;

  width = w;
  height = h;

  if ( h==0 )
    _d->m_aspect = ( GLdouble ) width;
  else
    _d->m_aspect = ( GLdouble ) width/ ( GLdouble ) height;


  glViewport ( 0,0,width,height );
  SetupProjection();

}

/** Setup the projection

   This method calls the SetupPerspectiveMethod and initializes the projection matrix
*/
void GLVisorBase::SetupProjection()
{
  glMatrixMode ( GL_PROJECTION );
  glLoadIdentity();
  SetupPerspective();
  glMatrixMode ( GL_MODELVIEW );
}


/** Setup the perspective by means of gluPerspective and gluLookAt functions*/
void GLVisorBase::SetupPerspective()
{
  gluPerspective ( 45.0,_d->m_aspect,1.0,60.0 );
  //glTranslatef ( _d->m_camera.x(),_d->m_camera.y(),_d->m_camera.z()-5.0 );
  float x=_d->m_camera.x();
  float y=_d->m_camera.y();
  float z=_d->m_camera.z();
  gluLookAt(x,y,z+5.0,x,y,z,0,1,0);
}

/** This method will call the RenderScene virtual method

   In most cases you should not override this method
*/
void GLVisorBase::paintGL()
{
  RenderScene();
}

/** \brief render the scene

   Reimplement in derived classes. Don't forget to call the base method
*/
void GLVisorBase::RenderScene ( GLenum mode )
{
  glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
}


/** Call the glMaterialfv function accordingly with the color associated to atom \a a*/
void GLVisorBase::GetMaterialColor ( const Atom &a )
{

  GLfloat color[4]={0.5,0.2,0.1,1.0};

  a.Color ( &color[0],&color[1],&color[2] );
  glMaterialfv ( GL_FRONT_AND_BACK, GL_AMBIENT, color );
  glMaterialfv ( GL_FRONT_AND_BACK, GL_DIFFUSE, color );
  return;



}

/** Call glColor3f accordingly to the color associated to the atom \a a*/
void GLVisorBase::GetColor ( const Atom& a )
{
  GLfloat color[3]={0.5,0.2,0.1};
  a.Color ( &color[0],&color[1],&color[2] );
  glColor3f ( color[0],color[1],color[2] );

}
/** manages mouse press events
  
  Don't forget to call the base method if you override this virtual method
  */
void GLVisorBase::mousePressEvent ( QMouseEvent* e )
{
  _d->m_oldpoint=e->pos();
  switch ( e->button() )
  {
    case Qt::LeftButton:
      _d->m_lgrabbed=true;
      _d->m_rgrabbed=false;
      _d->m_mgrabbed=false;
      selectControl ( e->pos() );
      break;
    case Qt::MidButton:
      _d->m_mgrabbed=true;
      _d->m_lgrabbed=false;
      _d->m_rgrabbed=false;
      break;
    case Qt::RightButton:
      _d->m_lgrabbed=false;
      _d->m_rgrabbed=true;
      _d->m_mgrabbed=false;
      break;
    default:
      break;
  }

}


void GLVisorBase::mouseReleaseEvent ( QMouseEvent* e )
{
  switch ( e->button() )
  {
    case Qt::LeftButton:
      _d->m_lgrabbed=false;
      break;
    case Qt::RightButton:
      _d->m_rgrabbed=false;
      break;
    case Qt::MidButton:
      _d->m_mgrabbed=false;
    default:
      break;
  }

  _d->m_bismoving=false;
  updateGL();

}

/** reimplement for custom handling of mouse move events*/
void GLVisorBase::mouseMoveEvent ( QMouseEvent* e )
{
  _d->m_oldpoint=e->pos();
  _d->m_bismoving=true;

  QGLWidget::mouseMoveEvent ( e );
}

/** reimplement for custom handling of mouse double click events*/
void GLVisorBase::mouseDoubleClickEvent ( QMouseEvent* e )
{
  QGLWidget::mouseDoubleClickEvent(e);
}

/** reimplement for custom handling of wheel events*/
void GLVisorBase::wheelEvent ( QWheelEvent* e )
{
  QGLWidget::wheelEvent(e);
}

/** show atom sequential numbers*/
void GLVisorBase::OnShowNumbers ( bool b )
{
  _d->m_bshownumbers=b;
  update();
}


/** call a QColorDialog widget and setup the background color*/
void GLVisorBase::OnChangeBackground()
{
  _d->m_bkcolor=QColorDialog::getColor ( _d->m_bkcolor );
  qglClearColor ( _d->m_bkcolor );
 update();
}

/** Show (or not) atom element symbols on screen*/
void GLVisorBase::OnShowSymbols( bool b )
{
  _d->m_bshowsymbols=b;
  update();
}

/** Show (or not) atom pdb info on screen*/
void GLVisorBase::OnShowPDBInfo( bool b )
{
  _d->m_bshowpdbinfo=b;
  update();
}


/** Override to do custom key press events handling*/
void GLVisorBase::keyPressEvent ( QKeyEvent* e )
{
  //let the parent manage also the event
  e->ignore();
}


/** \brief do initialization of the visor
    Override to do any necessary initalization
*/
void GLVisorBase::Initialize()
{

}

/** Select the currently active Molecule \a mol*/
void GLVisorBase::SelectMolecule ( size_t mol )
{
  _d->m_currentmolecule=mol;
  update();
}

/** Change the molecular rendering mode*/
void GLVisorBase::OnChangeGraphMode ( int mode )
{
  SetGraphMode ( static_cast<graphmode> ( mode ) );

}

void GLVisorBase::RenderText ( int x, int y, const QString& str )
{
  RenderText ( x,y,str,_d->m_glfont );
}

void GLVisorBase::RenderText ( int x, int y, const QString& str, const QFont& font )
{
  if ( _d->m_vectorgraphics )
  {
    gl2psText ( str.toStdString().c_str(),font.family().toStdString().c_str(),font.pointSize() );
  }
  else
    renderText ( x,y,str,font );
}

void GLVisorBase::RenderText ( double x, double y, double z, const QString& str )
{
  RenderText ( x,y,z,str,_d->m_glfont );
}

void GLVisorBase::RenderText ( double x, double y, double z, const QString& str, const QFont& font )
{
  if ( _d->m_vectorgraphics )
  {
    glRasterPos3d ( x,y,z );
    gl2psText ( str.toStdString().c_str(),font.family().toStdString().c_str(),font.pointSize() );

  }
  else
    renderText ( x,y,z,str,font );
}


/** Prepare the visor for exporting of the scene to a vector graphics file*/
void GLVisorBase::SetVectorGraphicsMode ( bool b )
{
  _d->m_vectorgraphics=b;
}


/** Use this function to ensure correct rendering of lines when exporting
vector pictures*/
void GLVisorBase::EnableLineStipple ( GLfloat width,GLint factor,GLushort pattern )
{
  if ( _d->m_vectorgraphics )
  {
    gl2psLineWidth ( width );
    glLineStipple ( factor,pattern );
    gl2psEnable ( GL_LINE_STIPPLE );

  }
  else
  {
    glLineWidth ( width );
    glEnable ( GL_LINE_STIPPLE );
    glLineStipple ( factor,pattern );

  }
  glDisable ( GL_LIGHTING );
}

/** Use this function to ensure correct rendering of lines when exporting
vector pictures*/
void GLVisorBase::DisableLineStipple()
{
  if ( _d->m_vectorgraphics )
  {
    gl2psDisable ( GL_LINE_STIPPLE );
  }
  else
  {
    glDisable ( GL_LINE_STIPPLE );
  }
  glEnable ( GL_LIGHTING );
}

/** \return a QPixmap of width \a w and height \a h with the
currently rendered scene*/
QPixmap GLVisorBase::Pixmap ( int w, int h )
{
  QPixmap p=renderPixmap ( w,h );
  resizeGL ( this->width(),this->height() );
  return p;

}

/** Set the number of spheres used in gluSphere and gluCylinder calls*/
void GLVisorBase::SetSpheres( int r )
{
  if ( r >= 1 )
    _d->m_spheres=r;
  else
    std::cerr << "not acceptable resolution";
}

/** \return the number of spheres used in gluSphere and gluCylinder calls*/
int GLVisorBase::Spheres() const
{
  return _d->m_spheres;
}

/** \brief Set the molecular rendering mode
   
  This method emits the graphmodeChanged signal
  and will update the screen
  */
void GLVisorBase::SetGraphMode ( graphmode mode )
{
  _d->m_graphmode=mode;
  emit graphmodeChanged(mode);
  update();
}

/** Set the color of the background*/
void GLVisorBase::SetBkColor ( const QColor& col )
{
  _d->m_bkcolor=col;
}

/** \return the font currently used for text rendering*/
const QFont& GLVisorBase::GLFont() const
{
  return _d->m_glfont;
}

/** \return the background color*/
QColor GLVisorBase::BkColor() const
{
  return _d->m_bkcolor;
}

/** Set the position of the camera*/
void GLVisorBase::SetCamera ( const Coordinate& c )
{
  _d->m_camera=c;
}
/** \return a const referece to a stl vector of trackball handlers*/
std::vector<MoleculeHandler>& GLVisorBase::Handlers()
{
  return _d->m_handlers;
}

/** \return a referece to a stl vector of trackball handlers*/
const std::vector<MoleculeHandler>& GLVisorBase::Handlers() const
{
  return _d->m_handlers;
}

/** \return true if mouse is currently moving*/
bool GLVisorBase::IsMouseMoving() const
{
  return _d->m_bismoving;
}

/** \return true if atom numbers are shown*/
bool GLVisorBase::ShowNumbers() const
{
  return _d->m_bshownumbers;
}

/** \return true if atom symbols are shown*/
bool GLVisorBase::ShowSymbols() const
{
  return _d->m_bshowsymbols;
}

/** \return true if atom pdb info is shown*/
bool GLVisorBase::ShowPDBInfo() const
{
  return _d->m_bshowpdbinfo;
}

/** \return true if dipole is shown*/
bool GLVisorBase::ShowDipole() const
{
  return _d->m_bshowdipole;
}

/** \return true if cell is shown*/
bool GLVisorBase::ShowCell() const
{
  return _d->m_bshowcell;
}

/** \return true if coordinate axis are shown*/
bool GLVisorBase::ShowAxis() const
{
  return _d->m_bshowaxis;
}

/** \return true if density is shown*/
bool GLVisorBase::ShowDensity() const
{
  return _d->m_bshowdensity;
}

/** \return the molecular rendering mode*/
GLVisorBase::graphmode GLVisorBase::GraphMode() const
{
  return _d->m_graphmode;
}

/** \return a Qt::MouseButtons combination of all grabbed mouse buttons*/
Qt::MouseButtons GLVisorBase::GrabbedButtons() const
{
  if ( _d->m_lgrabbed && _d->m_rgrabbed )
  {
     return ( Qt::LeftButton | Qt::RightButton );
  }
  
  if ( _d->m_lgrabbed )
  {
    return Qt::LeftButton;
  }
  
  if ( _d->m_rgrabbed )
  {
    return Qt::RightButton;
  }

  return ( Qt::LeftButton | Qt::RightButton );
  
}

/** \return the mouse position in widget Coordinates*/
const QPoint& GLVisorBase::MousePosition() const
{
  return _d->m_oldpoint;
}

/** \brief manage selection of  opengl rendered controls

 Reimplement this class if you want to use OpenGL rendered screen controls
 */
void GLVisorBase::selectControl ( const QPoint& pos )
{
}

QMenu* GLVisorBase::GetMenu() { return static_cast<QMenu*> ( m_optionsmenu ); }



