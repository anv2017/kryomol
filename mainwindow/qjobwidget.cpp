/*****************************************************************************************
                            qjobwidget.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "qjobwidget.h"

#include "world.h"

QJobWidget::QJobWidget( kryomol::World* world, QWidget* parent ) : QSplitter (parent), m_world (world)
{
}

QJobWidget::~QJobWidget()
{
}
