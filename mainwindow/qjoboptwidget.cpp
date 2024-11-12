/*****************************************************************************************
                            qjoboptwidget.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "qjoboptwidget.h"
#include "qconvwidget.h"

#include "world.h"
#include "glvisor.h"
#include "molecule.h"
#include "frame.h"
#include "kryovisor.h"

#include <QDockWidget>

QJobOptWidget::QJobOptWidget(QWidget* parent ) : QJobWidget (parent)
{
        m_world = new kryomol::World(this,kryomol::World::optvisor);
}

QJobOptWidget::~QJobOptWidget()
{
}


void QJobOptWidget::InitWidgets()
{
    World()->Visor()->Initialize();

    this->setCentralWidget(m_world->Visor());
    InitCommonWidgets();

    QDockWidget* cdock = new QDockWidget(this);
    m_convwidget = new QConvWidget (m_tabwidget);
    m_tabwidget->addTab(m_convwidget,"Opt");
    cdock->setWidget(m_tabwidget);
    cdock->setAllowedAreas(Qt::RightDockWidgetArea);
    this->addDockWidget(Qt::RightDockWidgetArea,cdock);

    // for(auto w : m_dockwidgets)
    // {
    //     this->tabifyDockWidget(cdock,w);
    // }

    // m_dockwidgets.insert(0,cdock);



    std::vector<kryomol::Frame> frames;
    std::vector<kryomol::Frame>::iterator mt;
    //remember to modify this horror
    for (mt=World()->Molecules().back().Frames().begin();mt!=World()->Molecules().back().Frames().end();mt++)
    {
        if ((mt->PotentialEnergy()) && (mt->RMSForce()))
            frames.push_back(*mt);
    }
    World()->Molecules().back().Frames().erase(World()->Molecules().back().Frames().begin(),World()->Molecules().back().Frames().end());
    for(mt=frames.begin();mt!=frames.end();mt++)
    {
        World()->Molecules().back().Frames().push_back(*mt);
    }

    m_convwidget->SetNData ( World()->Molecules().back().Frames().size() );
    double* energies=m_convwidget->GetEnergies();
    double* s2= m_convwidget->GetS2();
    double* rmsforce=m_convwidget->GetRMSForces();
    double* maximumforce=m_convwidget->GetMaximumForces();
    double* rmsdisplacement=m_convwidget->GetRMSDisplacements();
    double* maximumdisplacement=m_convwidget->GetMaximumDisplacements();
    m_convwidget->SetEnergyLevel ( World()->Molecules().back().GetEnergyLevel().c_str() );

    if ( energies )
    {
        std::vector<kryomol::Frame>::iterator mt;
        mt=World()->Molecules().back().Frames().begin();
        m_convwidget->SetThreshold ( mt->GetThreshold() );

        int i=0;
        mt=World()->Molecules().back().Frames().begin();
        for ( i=0;mt!=World()->Molecules().back().Frames().end();mt++,i++ )
        {
            energies[i]=mt->GetEnergy();
            s2[i]=mt->GetS2();
            rmsforce[i]=mt->GetRMSForce();
            maximumforce[i]=mt->GetMaximumForce();
            rmsdisplacement[i]=mt->GetRMSDisplacement();
            maximumdisplacement[i]=mt->GetMaximumDisplacement();
        }

    }

    //Initialize the visor and actions of the widget
    (static_cast<kryomol::KryoVisorOpt*> ( World()->Visor() ) )->setForceScale(sqrt(m_convwidget->GetForceScale()));
    m_convwidget->SetupCurves();
    connect ( m_convwidget,SIGNAL ( selectedPoint ( size_t) ),World(),SLOT ( SelectFrame(size_t ) ) );
    connect ( World(),SIGNAL ( currentFrame(size_t ) ),m_convwidget,SLOT ( OnSelectedPoint ( size_t ) ) );
    m_convwidget->OnSelectedPoint ( World()->CurrentMolecule()->CurrentFrameIndex());
    connect ( m_convwidget,SIGNAL ( forcescale ( float ) ),World()->Visor(),SLOT ( OnForceScale ( float ) ) );
    connect ( m_convwidget,SIGNAL ( showforces( bool ) ),World()->Visor(),SLOT ( OnShowForces( bool ) ) );

}
