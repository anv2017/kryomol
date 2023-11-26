/*****************************************************************************************
                            qcwtcombobox.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QCWTCOMBOBOX_H
#define QCWTCOMBOBOX_H

#include <QComboBox>

#include "export.h"

/**
  A combo box with custom what this*/
class KRYOMOL_API QCWTComboBox : public QComboBox
{
  public:
      QCWTComboBox(QWidget* parent = 0);
      virtual ~QCWTComboBox();
      bool customWhatsThis() const { return true; }
      protected:
      virtual void mousePressEvent(QMouseEvent* e);
  
  
};

#endif

