/*****************************************************************************************
                            nmrdevapp.h  -  description
                             -------------------
    copyright            : (C) 2011 by Armando Navarro-Vazquez and Noa Campos-Lopez
    email                : armando.navarro@uvigo.es, noa.campos@uvigo.es
******************************************************************************************/

#ifndef NMRDEVAPP_H
#define NMRDEVAPP_H
#include <QApplication>
#include <QDir>
#include <vector>
#include "export.h"

namespace nmrdev
{

class Plugin;
class PluginLoader;
class NMRDevApplicationPrivate;

/** \brief The main NMRDev application
*/
class NMRDEV_API NMRDevApplication : public QApplication
{
  public:
  NMRDevApplication ( int & argc, char ** argv);
  NMRDevApplication ( int & argc, char ** argv, bool GUIenabled);
  ~NMRDevApplication();
  
  static QString File();
  static void SetFile(const QString& file);
  static void SetFile(const char* file);
  
  std::vector<nmrdev::Plugin*>& Plugins();
  void LoadPlugins();
  void InitializePlugins();
  QString Authors() const;
  QString Version() const; 
  QString URL() const;
  protected:
  virtual bool notify(QObject* receiver, QEvent* e);
  
  private:
  NMRDevApplicationPrivate* _d;
};
}

#endif
