/*****************************************************************************************
                            scripter.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "scripter.h"
#include "qryomolscriptable.h"

#include <QLineEdit>
#include <QtScript>

using namespace kryomol;

Scripter::Scripter()
{
  m_commandline = new QLineEdit();
  connect(m_commandline,SIGNAL(returnPressed()),this,SLOT(OnEvaluateCommandLine()));
  m_engine = new QScriptEngine();
  m_scriptable = new KryoMolScriptable();
  QScriptValue worldproxy=m_engine->newQObject(m_scriptable);
  m_engine->globalObject().setProperty("world",worldproxy);
}

QLineEdit* Scripter::CommandLine()
{
  return m_commandline;
}

void Scripter::SetWorld(World* w)
{
  m_scriptable->SetWorld(w);
}

void Scripter::SetApplication(KryoMolApplication* app)
{
  m_scriptable->SetApplication(app);
}
void Scripter::OnEvaluateCommandLine()
{
  QString text=m_commandline->text();
  m_engine->evaluate(text);
}
