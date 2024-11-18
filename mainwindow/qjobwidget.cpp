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
#include "glvisor.h"

#include "qpdbcontrol.h"
#include "qmeasurewidget.h"
#include "qorbitalwidget.h"
#include "qmolecularlistcontrol.h"

#include <QDockWidget>

QJobWidget::QJobWidget(QWidget* parent) : QMainWindow (parent), m_world (nullptr)
{
    m_tabwidget = new QTabWidget();
    //This is necessary to the right docked widget occupy the whole height
    this->setCorner ( Qt::BottomRightCorner, Qt::RightDockWidgetArea );
}

QJobWidget::~QJobWidget()
{
    delete m_world;
}


void QJobWidget::showEvent(QShowEvent *event)
{
    for(auto w : m_dockwidgets)
    {
        w->show();
    }
}

void QJobWidget::hideEvent(QHideEvent *event)
{
    for(auto w : m_dockwidgets)
    {
        w->hide();
    }

}

void QJobWidget::InitCommonWidgets()
{
    /*QDockWidget* pdbdock = new QDockWidget(this);
    pdbdock->setWidget( new QPDBControl(pdbdock));
    m_dockwidgets.push_back(pdbdock);

    QDockWidget* mdock = new QDockWidget(this);
    pdbdock->setWidget( new QMeasureWidget(mdock));
    m_dockwidgets.push_back(mdock);

    QDockWidget* orbdock = new QDockWidget(this);
    orbdock->setWidget( new QOrbitalWidget(this));
    m_dockwidgets.push_back(orbdock);

    QDockWidget* mddock= new QDockWidget(this);
    mddock->setWidget( new QMolecularListControl(mddock));
    m_dockwidgets.push_back(mddock);


    for(auto w : m_dockwidgets)
    {
        w->setAllowedAreas(Qt::RightDockWidgetArea);
        this->addDockWidget(Qt::RightDockWidgetArea,w);
    }*/
    QMolecularListControl* conf= new QMolecularListControl(m_tabwidget);
    conf->SetWorld(m_world);
    conf->Init();
    m_tabwidget->addTab(conf,"Conformers");

    m_tabwidget->addTab(new QPDBControl(m_tabwidget),"PDB Info");
    QMeasureWidget* mw=new QMeasureWidget(m_tabwidget);
    m_tabwidget->addTab(mw,"Measure");
    m_tabwidget->addTab(new QOrbitalWidget(m_tabwidget),"Density");
    //Add measure connections
    //CONNECTIONS
    connect( m_world->Visor(), SIGNAL ( distance ( QString& ) ), mw, SLOT ( OnWriteDistance ( QString& ) ) );
    //connect( m_world->Visor(), SIGNAL ( angle ( QString& ) ), m_measures, SLOT ( OnWriteAngle ( QString& ) ) );
    //connect( m_world->Visor(), SIGNAL ( dihedral ( QString& ) ), m_measures, SLOT ( OnWriteDihedral ( QString& ) ) );
    //connect(m_world, SIGNAL ( currentFrame ( size_t ) ), this, SLOT ( OnUpdateMeasures ( size_t )) );*/
    //connect(m_measures,SIGNAL(clearAll()),m_world->Visor(), SLOT(OnClearMeasures()));
    //connect(m_measures,SIGNAL(distanceChange(int)),m_world->Visor(),SLOT(OnDistanceChange(int)));
    //connect(m_measures,SIGNAL(angleChange(int)),m_world->Visor(),SLOT(OnAngleChange(int)));
    //connect(m_measures,SIGNAL(dihedralChange(int)),m_world->Visor(),SLOT(OnDihedralChange(int)));
    //connect(m_measures, SIGNAL(showDistances(bool)),m_world->Visor(),SLOT(OnShowDistances(bool)));


}
