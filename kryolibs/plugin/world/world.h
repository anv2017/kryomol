/*****************************************************************************************
                            world.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef WORLD_H
#define WORLD_H

#include "thermo.h"
#include "export.h"
//Added by qt3to4:
#include <QMouseEvent>
#include <QDropEvent>
#include <QDragEnterEvent>
#include <QObject>

class QGLWidget;

namespace kryomol
{
  class Plugin;
  class GLVisor;
  class WorldPrivate;
  class GLVisorBase;
  class KryoVisor;
  
  /**
  \brief simulation world of KryoMol

  The World class represents the  world simulation
  containing all loaded molecules and thermodynamic parameters
  */
  class KRYOMOL_API World : public QObject, public Thermostat
  {
      Q_OBJECT
    public:
      enum VisorType {glvisor, optvisor, freqvisor, opticalvisor, kryovisor};
      World ( bool bGUI=false );
      World ( QWidget* parent, VisorType vtype, const QGLWidget* shareWidget=0, Qt::WindowFlags f=0 );
      ~World();
      virtual void Initialize();
      void SetCurrentPlugin ( Plugin* p );
      Plugin* CurrentPlugin();
      const Plugin* CurrentPlugin() const;
      virtual void SetTemperature ( double t );
      virtual void SetEnsamble ( ensamble e );
      GLVisor* Visor();
      const GLVisor* Visor() const;
      void Clear();
      std::vector<Molecule>& Molecules();
      const std::vector<Molecule>& Molecules() const;
      const Molecule* CurrentMolecule() const;
      Molecule* CurrentMolecule();
      size_t CurrentMoleculeIndex() const;

    signals:
      /** emitted when user changes temperatue or the kind of thermodynamic ensamble*/
      void thermostatChanged();
      /** emited when the current frame changes*/
      void currentFrame ( size_t );
      /** emitted when changed current molecule */
      void currentMolecule(size_t );

    public slots:
      void SelectFrame ( size_t );
      void SelectMolecule ( size_t );

    private:
        #ifdef __GNUC__
        #warning implement d pointer
        #endif
      Plugin* m_currentplugin;
      GLVisor* m_visor;
      std::vector<Molecule> m_molecules;
      size_t m_currentmolecule;

  };
}
#endif
