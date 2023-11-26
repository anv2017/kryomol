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

QJobOptWidget::QJobOptWidget( kryomol::World* world, QWidget* parent ) : QJobWidget (world, parent)
{
    Init();
}

QJobOptWidget::~QJobOptWidget()
{
}


void QJobOptWidget::Init()
{
    kryomol::World* world = GetWorld();

    m_convwidget = new QConvWidget (this);

    std::vector<kryomol::Frame> frames;
    std::vector<kryomol::Frame>::iterator mt;
    //remember to modify this horror
    for (mt=world->Molecules().back().Frames().begin();mt!=world->Molecules().back().Frames().end();mt++)
    {
        if ((mt->PotentialEnergy()) && (mt->RMSForce()))
            frames.push_back(*mt);
    }
    world->Molecules().back().Frames().erase(world->Molecules().back().Frames().begin(),world->Molecules().back().Frames().end());
    for(mt=frames.begin();mt!=frames.end();mt++)
    {
        world->Molecules().back().Frames().push_back(*mt);
    }

    m_convwidget->SetNData ( world->Molecules().back().Frames().size() );
    double* energies=m_convwidget->GetEnergies();
    double* s2= m_convwidget->GetS2();
    double* rmsforce=m_convwidget->GetRMSForces();
    double* maximumforce=m_convwidget->GetMaximumForces();
    double* rmsdisplacement=m_convwidget->GetRMSDisplacements();
    double* maximumdisplacement=m_convwidget->GetMaximumDisplacements();
    m_convwidget->SetEnergyLevel ( world->Molecules().back().GetEnergyLevel().c_str() );

    if ( energies )
    {
        std::vector<kryomol::Frame>::iterator mt;
        mt=world->Molecules().back().Frames().begin();
        m_convwidget->SetThreshold ( mt->GetThreshold() );

        int i=0;
        mt=world->Molecules().back().Frames().begin();
        for ( i=0;mt!=world->Molecules().back().Frames().end();mt++,i++ )
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
    world->Visor()->Initialize();
    (static_cast<kryomol::KryoVisorOpt*> ( world->Visor() ) )->setForceScale(sqrt(m_convwidget->GetForceScale()));
    m_convwidget->SetupCurves();
    connect ( m_convwidget,SIGNAL ( selectedPoint ( size_t) ),world,SLOT ( SelectFrame(size_t ) ) );
    connect ( world,SIGNAL ( currentFrame(size_t ) ),m_convwidget,SLOT ( OnSelectedPoint ( size_t ) ) );
    m_convwidget->OnSelectedPoint ( world->CurrentMolecule()->CurrentFrameIndex());
    connect ( m_convwidget,SIGNAL ( forcescale ( float ) ),world->Visor(),SLOT ( OnForceScale ( float ) ) );
    connect ( m_convwidget,SIGNAL ( showforces( bool ) ),world->Visor(),SLOT ( OnShowForces( bool ) ) );

    //Add the visor and the ConvWidget to the splitter
    this->addWidget(world->Visor());
    this->addWidget(m_convwidget);

    SetWorld(world);
}
