/*****************************************************************************************
                            qjobspwidget.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "qjobspwidget.h"
#include "world.h"
#include "kryovisor.h"

QJobSpWidget::QJobSpWidget(QWidget* parent ) : QJobWidget (parent)
{
    m_world = new kryomol::World(this,kryomol::World::glvisor);
}

QJobSpWidget::~QJobSpWidget()
{
}


void QJobSpWidget::InitWidgets()
{
    this->World()->Initialize();
    this->setCentralWidget(m_world->Visor());
    InitCommonWidgets();
}


