/*****************************************************************************************
                            glvisorbase.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef GLVISORBASE_H
#define GLVISORBASE_H

#include <QtOpenGL/QGLWidget>
#include <QObject>
//Added by qt3to4:
#include <QWheelEvent>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QMenu>
#include <qwidget.h>
#include <qcolor.h>
#include <qtimeropenmenu.h>

#include "molecule.h"
#include "export.h"
#include "handler.h"
#include "kryomolcore_export.h"


struct cartesian {GLdouble x; GLdouble y; GLdouble z;};

class QInfoDialog;
class QMenu;

namespace kryomol
{
   
   class GLVisorBasePrivate;

  /** @brief base class for OpenGL rendered widget

  This QGLWidget derived class do basic rendering of molecules in several modes (Wireframe, Sticks and CPK)
  It allows exporting of vector pictures through the use of the gl2ps library.
  To ensure proper vector picture exporting use the RenderText, EnableLineStipple and DisableLineStipple methods
  */
  class KRYOMOL_API GLVisorBase : public QGLWidget
  {
     public:
      /** the molecular rendering modes*/
      enum graphmode { WIREFRAME=0, STICKS, CPK };

      Q_OBJECT

    public:
      GLVisorBase ( QWidget* parent=0, const QGLWidget* shareWidget=0,Qt::WindowFlags f=0 );
      virtual ~GLVisorBase();
      virtual void Initialize();
      void SetGraphMode ( graphmode mode );
      graphmode GraphMode() const;
      QColor BkColor() const ;
      void SetBkColor ( const QColor& bk );
      /** Render text, in screen coordinates, either on screen or vector picture*/
      void RenderText ( int x, int y, const QString& str, const QFont& );
      /** Render text, in screen coordinates, either on screen or vector picture using widget font*/
      void RenderText ( int x, int y, const QString& str );
      /** Render text, in model coordinates, either on screen or vector picture using widget font*/
      void RenderText ( double x, double y, double z, const QString& str,const QFont& );
      /** Render text, in model coordinates, either on screen or vector picture using widget font*/
      void RenderText ( double x, double y, double z, const QString& str );
      void EnableLineStipple ( GLfloat width=1,GLint factor=1,GLushort pattern=0xCCCC );
      void DisableLineStipple();
      QPixmap Pixmap ( int w, int h );
      void SetSpheres ( int r );
      int Resolution() const;
      void SetCamera ( const Coordinate& c );
      const QFont& GLFont() const;
      Qt::MouseButtons GrabbedButtons() const;
      bool ShowNumbers() const;
      bool ShowSymbols() const;
      bool ShowPDBInfo() const;
      bool ShowDipole() const;
      bool ShowCell() const;
      bool ShowAxis() const;
      bool ShowDensity() const;

    public slots:
      void OnShowNumbers ( bool b );
      void OnChangeBackground();
      void OnShowSymbols ( bool b );
      void OnShowPDBInfo(bool b);
      void SelectMolecule ( size_t current );
      int Spheres() const;

    protected:
      void GetColor ( const Atom &a );
      void GetMaterialColor ( const Atom& a );
      void SetVectorGraphicsMode ( bool b );
      void SetupPerspective();
      void SetupProjection();
      const std::vector<MoleculeHandler>& Handlers() const;
      std::vector<MoleculeHandler>& Handlers();
      bool IsMouseMoving() const;
      const QPoint& MousePosition() const;
      QMenu* GetMenu();

    protected:
      virtual void selectControl ( const QPoint& pos );
      virtual void initializeGL();
      virtual void resizeGL ( int w, int h );
      virtual void paintGL();
      virtual void wheelEvent ( QWheelEvent* e );
      virtual void RenderScene ( GLenum mode=GL_RENDER );
      virtual void OnSelection ( QMouseEvent* e ) {}
      virtual void mouseMoveEvent ( QMouseEvent * e );
      virtual void mousePressEvent ( QMouseEvent* e );
      virtual void mouseReleaseEvent ( QMouseEvent* e );
      virtual void mouseDoubleClickEvent ( QMouseEvent* e );
      virtual void keyPressEvent ( QKeyEvent* e );

    private slots:
      void OnChangeGraphMode ( int mode );

    signals:
      /** emit a message to be displayed, for instance in the status bar*/
      void  text ( const QString& );
      /** emitted when the molecular rendering mode has been changed*/
      void graphmodeChanged ( kryomol::GLVisorBase::graphmode );

    private:
      void Init();

    private:
      QWidget* m_parent;
      GLfloat m_range;
      GLVisorBasePrivate* _d;
      QTimerOpenMenu* m_optionsmenu;

  };

}

#endif
