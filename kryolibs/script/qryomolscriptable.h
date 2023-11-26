/*****************************************************************************************
                            qryomolscriptable.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/


#ifndef QRYOMOLSCRIPTABLE_H
#define QRYOMOLSCRIPTABLE_H

#include <QtScript>
#include "scriptexport.h"

namespace kryomol
{
class World;
class KryoMolApplication;

class SCRIPT_API  KryoMolScriptable : public QObject, protected QScriptable
{
  Q_OBJECT
  public:
  KryoMolScriptable();
  void SetWorld(World* w);
  void SetApplication(KryoMolApplication* app);
  public slots:
  void popmessage(QString text);
  void superimpose(QString text);
  void eckarttransform(QString text);
  private:
  World* m_world;
  KryoMolApplication* m_app;
};
}
#endif



