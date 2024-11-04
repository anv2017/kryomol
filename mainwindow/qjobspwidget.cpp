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
    World()->Visor()->Initialize();
    InitCommonWidgets();

    for(auto it=m_dockwidgets.begin()+1;it!=m_dockwidgets.end();++it)
    {
        this->tabifyDockWidget(m_dockwidgets.front(),*it);
    }


}


