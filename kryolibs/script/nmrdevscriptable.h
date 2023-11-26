/*****************************************************************************************
                            nmrdevscriptable.h  -  description
                             -------------------
    copyright            : (C) 2011 by Armando Navarro-Vazquez and Noa Campos-Lopez
    email                : armando.navarro@uvigo.es, noa.campos@uvigo.es
******************************************************************************************/


#ifndef NMRDEVSCRIPTABLE_H
#define NMRDEVSCRIPTABLE_H

#include <QtScript>
#include "scriptexport.h"

namespace nmrdev
{
class World;
class NMRDevApplication;
class SCRIPT_API  NMRDevScriptable : public QObject, protected QScriptable
{
  Q_OBJECT
  public:
  NMRDevScriptable();
  void SetWorld(World* w);
  void SetApplication(NMRDevApplication* app);
  public slots:
  void popmessage(QString text);
  void superimpose(QString text);
  void eckarttransform(QString text);
  private:
  World* m_world;
  NMRDevApplication* m_app;
};
}
#endif



