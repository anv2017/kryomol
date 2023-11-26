/*****************************************************************************************
                            qryomolapp.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QRYOMOLAPP_H
#define QRYOMOLAPP_H
#include <QApplication>
#include <QDir>
#include <vector>
#include "export.h"

namespace kryomol
{

class Plugin;
class PluginLoader;
class KryoMolApplicationPrivate;

/** \brief The main QryoMol application
*/
class KRYOMOL_API KryoMolApplication : public QApplication
{
  public:
  KryoMolApplication ( int & argc, char ** argv);
  KryoMolApplication ( int & argc, char ** argv, bool GUIenabled);
  ~KryoMolApplication();
  
  static QString File();
  static void SetFile(const QString& file);
  static void SetFile(const char* file);
  
  std::vector<kryomol::Plugin*>& Plugins();
  void LoadPlugins();
  void InitializePlugins();
  QString Authors() const;
  QString Version() const; 
  QString URL() const;
  protected:
  virtual bool notify(QObject* receiver, QEvent* e);
  
  private:
  KryoMolApplicationPrivate* _d;
};
}

#endif
