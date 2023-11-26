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

QJobSpWidget::QJobSpWidget(kryomol::World* world, QWidget* parent ) : QJobWidget (world, parent)
{
    Init();
}

QJobSpWidget::~QJobSpWidget()
{
}


void QJobSpWidget::Init()
{
    kryomol::World* world = GetWorld();

    world->Visor()->Initialize();

    //Add the visor and the FreqWidget to the splitter
    this->addWidget(world->Visor());

    SetWorld(world);
}


