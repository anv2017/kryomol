/***************************************************************************
                          qtimeropenmenu.h  -  description
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

#ifndef QTIMEROPENMENU_H
#define QTIMEROPENMENU_H

#include <QMenu>
//Added by qt3to4:
#include <QEvent>

/**
@author Armando Navarro-Vazquez
*/

class QTimer;
class QTimerOpenMenu : public QMenu
{
  Q_OBJECT
public:
    QTimerOpenMenu(QWidget* parent=0);

    ~QTimerOpenMenu();
protected:
virtual bool eventFilter(QObject*,QEvent* );
public slots:
/** execute exec(QCursor::pos)*/
  void start(int time= 5 /*ms*/);   
private slots:
  void open();
private:
  QTimer* m_timer;
  QPoint m_pos;
  

};

#endif
