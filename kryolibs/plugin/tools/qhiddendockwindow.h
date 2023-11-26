/*****************************************************************************************
                            qhiddendockwindow.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QHIDDENDOCKWINDOW_H
#define QHIDDENDOCKWINDOW_H


#include <QDockWidget>
//Added by qt3to4:
#include <QHideEvent>
#include "export.h"


/**
a specialized class for hidden aware dock windows
*/


class KRYOMOL_API QHiddenDockWidget : public QDockWidget
{
Q_OBJECT

public:
  QHiddenDockWidget(const QString& title, QWidget* parent =0, Qt::WindowFlags flags = 0);
  QHiddenDockWidget (QWidget * parent = 0, Qt::WindowFlags f = 0 );
  ~QHiddenDockWidget();

signals:
  void hidden();

protected:
  virtual void hideEvent(QHideEvent* );

};

#endif
