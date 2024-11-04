#include "qjobdynwidget.h"
#include "world.h"
#include "glvisor.h"
#include "molecule.h"
#include "qconvwidget.h"
#include "kryovisor.h"

#include <QDockWidget>

QJobDynWidget::QJobDynWidget(QWidget *parent) : QJobWidget(parent)
{
    m_world = new kryomol::World(kryomol::World::glvisor);
}

void QJobDynWidget::InitWidgets()
{
    InitCommonWidgets();
    this->setCentralWidget(m_world->Visor());

    QDockWidget* dyndock = new QDockWidget(this);
    m_dynwidget = new QConvWidget (dyndock,false,false);
    dyndock->setWidget(m_dynwidget);
    dyndock->setAllowedAreas(Qt::RightDockWidgetArea);

    for(auto w : m_dockwidgets)
    {
        this->tabifyDockWidget(dyndock,w);
    }

    m_dockwidgets.insert(0,dyndock);

    //delete last frame if we do not have potential energy
    if ( !World()->Molecules().back().Frames().back().PotentialEnergy() )
    {
        World()->Molecules().back().Frames().pop_back();
    }


    m_dynwidget->SetNData ( World()->Molecules().back().Frames().size() );
    double* energies=m_dynwidget->GetEnergies();
    m_dynwidget->SetEnergyLevel ( World()->Molecules().back().GetEnergyLevel().c_str() );

    if ( energies )
    {
        std::vector<kryomol::Frame>::iterator mt;
        mt=World()->Molecules().back().Frames().begin();
        m_dynwidget->SetThreshold ( mt->GetThreshold() );

        int i=0;
        mt=World()->Molecules().back().Frames().begin();
        for ( i=0;mt!=World()->Molecules().back().Frames().end();mt++,i++ )
        {
            energies[i]=mt->GetEnergy();
        }

    }

    //Initialize the visor and actions of the widget
    World()->Visor()->Initialize();
    //(static_cast<kryomol::KryoVisorOpt*> ( world->Visor() ) )->setForceScale(sqrt(m_convwidget->GetForceScale()));
    m_dynwidget->SetupCurves();
    connect ( m_dynwidget,SIGNAL ( selectedPoint ( size_t) ),World(),SLOT ( SelectFrame(size_t ) ) );
    connect ( World(),SIGNAL ( currentFrame(size_t ) ),m_dynwidget,SLOT ( OnSelectedPoint ( size_t ) ) );
    m_dynwidget->OnSelectedPoint ( World()->CurrentMolecule()->CurrentFrameIndex());
}
