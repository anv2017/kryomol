/*****************************************************************************************
                            qtimeropenmenu.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <iostream>

#include <qtimer.h>
#include <qcursor.h>
#include "qtimeropenmenu.h"
//Added by qt3to4:
#include <QMenu>
#include <QEvent>

QTimerOpenMenu::QTimerOpenMenu(QString& name, QWidget* parent)
    : QMenu(name, parent)
{
  m_timer=new QTimer(this);
  parent->installEventFilter(this);
  connect(m_timer,SIGNAL(timeout()),this,SLOT(open()));

}


QTimerOpenMenu::~QTimerOpenMenu()
{}

void QTimerOpenMenu::start(int time /*=5*/)
{
  m_pos=QCursor::pos();
  m_timer->start(time);

}

void QTimerOpenMenu::open()
{
  this->exec(QCursor::pos());
}

/* calling mouseMoveEvent(QMouseEvent*) does not work till the menu is shown*/
/** This introduces a lot of overhead so try to search for another solution*/
bool QTimerOpenMenu::eventFilter(QObject* watched,QEvent* e)
{

  if( e->type()==QEvent::MouseMove && m_timer->isActive() )
  {
    m_timer->stop();
  }
  return false;
}


