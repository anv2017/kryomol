#include "qjobdynwidget.h"
#include "world.h"
#include "glvisor.h"
#include "molecule.h"
#include "qconvwidget.h"
#include "kryovisor.h"

QJobDynWidget::QJobDynWidget( kryomol::World* world, QWidget *parent) : QJobWidget(world,parent)
{
  Init();
}

void QJobDynWidget::Init()
{
    kryomol::World* world = GetWorld();

    m_dynwidget = new QConvWidget (this,false,false);

    //delete last frame if we do not have potential energy
    if ( ! world->Molecules().back().Frames().back().PotentialEnergy() )
    {
       world->Molecules().back().Frames().pop_back();
    }


    m_dynwidget->SetNData ( world->Molecules().back().Frames().size() );
    double* energies=m_dynwidget->GetEnergies();
    m_dynwidget->SetEnergyLevel ( world->Molecules().back().GetEnergyLevel().c_str() );

    if ( energies )
    {
        std::vector<kryomol::Frame>::iterator mt;
        mt=world->Molecules().back().Frames().begin();
        m_dynwidget->SetThreshold ( mt->GetThreshold() );

        int i=0;
        mt=world->Molecules().back().Frames().begin();
        for ( i=0;mt!=world->Molecules().back().Frames().end();mt++,i++ )
        {
            energies[i]=mt->GetEnergy();
        }

    }

    //Initialize the visor and actions of the widget
    world->Visor()->Initialize();
    //(static_cast<kryomol::KryoVisorOpt*> ( world->Visor() ) )->setForceScale(sqrt(m_convwidget->GetForceScale()));
    m_dynwidget->SetupCurves();
    connect ( m_dynwidget,SIGNAL ( selectedPoint ( size_t) ),world,SLOT ( SelectFrame(size_t ) ) );
    connect ( world,SIGNAL ( currentFrame(size_t ) ),m_dynwidget,SLOT ( OnSelectedPoint ( size_t ) ) );
    m_dynwidget->OnSelectedPoint ( world->CurrentMolecule()->CurrentFrameIndex());

    //Add the visor and the ConvWidget to the splitter
    this->addWidget(world->Visor());
    this->addWidget(m_dynwidget);

    SetWorld(world);
}
