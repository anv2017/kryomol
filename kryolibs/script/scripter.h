/*****************************************************************************************
                            scripter.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef SCRIPTER_H
#define SCRIPTER_H

#include <QObject>

#include "scriptexport.h"

class QWidget;
class QLineEdit;
class QScriptEngine;

namespace kryomol
{
class World;
class KryoMolScriptable;
class KryoMolApplication;

class SCRIPT_API Scripter : public QObject
{
    Q_OBJECT
public:
    Scripter();
    QLineEdit* CommandLine();
    void SetWorld(World* w);
    void SetApplication(KryoMolApplication* app);
private slots:
    void OnEvaluateCommandLine();
private:
    QLineEdit* m_commandline;
    QScriptEngine* m_engine;
    KryoMolScriptable* m_scriptable;
};
}
#endif
