/*****************************************************************************************
                            qtimeropenmenu.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QTIMEROPENMENU_H
#define QTIMEROPENMENU_H

#include <qmenu.h>
//Added by qt3to4:
#include <QEvent>

class QTimer;
class  QTimerOpenMenu : public QMenu
{
  Q_OBJECT
public:
    QTimerOpenMenu(QString& name, QWidget* parent=0);

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
