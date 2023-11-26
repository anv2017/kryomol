/*****************************************************************************************
                            qhiddendockwindow.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "qhiddendockwindow.h"
//Added by qt3to4:
#include <QHideEvent>


QHiddenDockWidget::QHiddenDockWidget ( const QString& title, QWidget * parent, Qt::WindowFlags f) : QDockWidget(title,parent,f)
{}

QHiddenDockWidget::QHiddenDockWidget(QWidget * parent, Qt::WindowFlags f) : QDockWidget(parent,f)
{}


QHiddenDockWidget::~QHiddenDockWidget()
{}

void QHiddenDockWidget::hideEvent(QHideEvent* e)
{
  emit hidden();
  QDockWidget::hideEvent(e);

}

