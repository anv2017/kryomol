/*****************************************************************************************
                            animation.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef ANIMATION_H
#define ANIMATION_H

#include <QObject>
#include <QTimer>

#include "molecule.h"
#include "coordinate.h"

class QTimer;

namespace kryomol {

    class Animation : public QObject
    {
      Q_OBJECT
    public:
      Animation (Molecule& molecule);
      ~Animation();
      void SetMode(int mode,int frame);
      int GetMode() const { return m_animationmode; }
      void SetFrame(int frame);
      void Start();
      void Stop();
      void Clear() { m_backup.clear(); }
      void Break();
      void BackupMolecule();
      void RestoreMolecule();
      bool IsActive() const;

    signals:
      void shot();
    private slots:
      void  OnAnimationTimer();
    private:
      Molecule& m_molecule;
      int m_animationmode;
      QTimer* m_timer;
      std::vector<Coordinate> m_backup;
      int m_interval;
      bool m_bfirst;
      int m_counter;
      int m_direction;
      /** store here the values for the distortions*/
      std::vector<int> m_previousframe;

    };
}
#endif
