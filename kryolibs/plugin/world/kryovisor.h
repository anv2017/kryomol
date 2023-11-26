/*****************************************************************************************
                            kryovisor.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef KRYOVISOR_H
#define KRYOVISOR_H

#include "glvisor.h"
//Added by qt3to4:
#include <QWheelEvent>
#include <QMouseEvent>
#include <QMenu>
#include <iostream>
#include "kryomolcore_export.h"

/** Spetialization for the KryoMol Plugin*/

class QMenu;

namespace kryomol
{
    class World;

    enum eSelectionMode { eNone,  eJ, eDistance, eAngle, eDihedral, eAnglePlanes, eRotateBond, eCreateBond,eDeleteBond, eAddGroup, eRenumberAll, eRenumberHeavy, eSetDihedral, eGromacs3dout };

    class Animation;

    class KRYOMOLCORE_EXPORT KryoVisor : public GLVisor
    {
     Q_OBJECT
    public:
       KryoVisor(World* world,QWidget* parent =0,const QGLWidget* shareWidget=0, Qt::WindowFlags f=0 );

       virtual ~KryoVisor();

      virtual void Initialize();
      void Center();

    signals:
      void selectedPoint(size_t );
      void playing(bool );

    public slots:
      void OnSelectionNone() { m_selectionmode=NONE; OnResetSelection(true);}
      void OnClearDistances();
      void OnMeasureDistances(bool b);
      void OnMeasureAngles(bool b);
      void OnMeasureDihedrals(bool b);
      void OnRotateBond(bool b);
      void OnSelectPoint(size_t current);
      void OnFirstFrame();
      void OnLastFrame();
      void OnPreviousFrame();
      void OnNextFrame();

    protected:
      virtual void OnSelection(QMouseEvent* e);
      virtual void RenderScreenText();
      virtual void RenderScene(GLenum mode);
      virtual void mousePressEvent(QMouseEvent* e);
      virtual void wheelEvent(QWheelEvent* e);
      virtual void mouseMoveEvent(QMouseEvent* e);
      virtual void RenderMolecule(size_t index,GLenum mode=GL_RENDER);
      void SetSelectionMode(selectionmode mode);
      selectionmode SelectionMode() const { return m_selectionmode; }

    private:
      /** Change hydrogen by a methyl group, this should be substitued for a more general function*/
      void ChangeHByMethyl(size_t i);
      void RenumberAtoms();
      void Gromacs3dout(size_t i, size_t j, size_t k, size_t l);
      std::vector<Coordinate> CalculateCellCube();

    protected slots:
      void OnResetSelection(bool repaint =true);
      void OnRecalculateConnectivity();
      void OnMoveToCentroid();
      void OnShowESPCharges();

    protected:
      void RefreshDistances();

    protected:
      World* m_world;
      std::vector<Molecule::pair> m_rpairs;
      std::vector<Molecule::triad> m_angletriads;
      std::vector<Molecule::quad> m_dihedralquads;
      QMenu* m_exportmenu;
      /** Show dipole on the screen*/
      bool m_bshowdipole;
      bool m_bshowcell;
      bool m_bshowaxis;
      /** Show bond lenghts on screen*/
      selectionmode m_selectionmode;
      QTimer* m_timer;
      int m_framespeed;
      std::vector< std::vector<Molecule::Coupling>::const_iterator > m_jpairs;
      unsigned int  m_nheavyatoms;
      std::vector<QRect> m_screencontrols;

    private:
      QMenu* m_distancemenu;
      QMenu* m_chargemenu;
      bool m_bshowespcharges;
      std::vector<size_t> m_selatoms;
      class KryoVisorPrivate;
      KryoVisorPrivate* _d;
    };


    class KRYOMOLCORE_EXPORT KryoVisorOpt : public KryoVisor
    {
        Q_OBJECT

    public:
      KryoVisorOpt(World* world,QWidget* parent=0, const QGLWidget* shareWidget=0, Qt::WindowFlags f=0 ) ;
      virtual ~KryoVisorOpt() {}
      virtual void Initialize();
      void setForceScale(int v) {m_forcescale=v;}

    public slots:
      void OnForceScale(float f) { m_forcescale=f; update(); }      
      void OnShowForces(bool b);
      void OnStartAnimation();
      void OnStopAnimation();

    protected:
      virtual void RenderScreenText();
      virtual void RenderScene(GLenum mode);
      virtual void RenderMolecule(size_t index, GLenum mode=GL_RENDER);
      virtual void wheelEvent(QWheelEvent* e);

    protected slots:
      void OnChangeFrame();

    protected:
      float m_forcescale;
      float m_forcemagnitude;
      bool m_bshowforces;

    };

    class QPNGMovie;
    class KRYOMOLCORE_EXPORT KryoVisorFreq : public KryoVisor
    {
      Q_OBJECT

    public:
      KryoVisorFreq(World* world,QWidget* parent=0, const QGLWidget* shareWidget=0, Qt::WindowFlags f=0 );
      virtual ~KryoVisorFreq();
      virtual void Initialize();
      bool IsPlaying() const;

    public slots:
      void SetMode(int mode,int frame);
      void OnDistortion(int);
      void OnStartAnimation();
      void OnStopAnimation();

    signals:
      /** emited when vibrational modes are animated or stopped*/
      void playing(bool);

    private slots:
      void OnPushImage();

    protected:
      virtual void RenderScreenText();      
      virtual void RenderScene(GLenum mode);
      virtual void RenderMolecule(size_t index,GLenum mode=GL_RENDER);
      virtual void wheelEvent(QWheelEvent* e);
      Animation* m_animation;

    };

    template <class X> class KryoVisorFactory
    {
      public:
      static X* BuildVisor(World* world,QWidget* parent)
      {
       return new X(parent);
      }
    };
}
#endif
