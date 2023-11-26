/*****************************************************************************************
                            qcwtcombobox.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "qcwtcombobox.h"
#include <QMouseEvent>
#include <QWhatsThis>

QCWTComboBox::QCWTComboBox(QWidget* parent) : QComboBox(parent)
{}

QCWTComboBox::~QCWTComboBox()
{}



void QCWTComboBox::mousePressEvent(QMouseEvent* e)
{
  if ( e->button() ==Qt::LeftButton)
  {
    if (QWhatsThis::inWhatsThisMode())
    {
      QVariant hdata=itemData(currentIndex());
      if ( hdata.isValid() )
        QWhatsThis::showText(QCursor::pos(),hdata.toString());
    }
    else
      QComboBox::mousePressEvent(e);
  }
  else
    QComboBox::mousePressEvent(e);
}
