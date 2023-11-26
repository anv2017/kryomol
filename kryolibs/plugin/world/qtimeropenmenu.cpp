/***************************************************************************
                          qtimeropenmenu.cpp  -  description
                             -------------------
    copyright            : (C) 2007 by A. Navarro-Vazquez
    email                : qoajnv@usc.es
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/


#include <iostream>

#include <qtimer.h>
#include <qcursor.h>
#include "qtimeropenmenu.h"

QTimerOpenMenu::QTimerOpenMenu(QWidget* parent)
    : QMenu(parent)
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
  m_timer->start(time,true);

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
    m_timer->stop();
  return false;
}


