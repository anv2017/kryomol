/*****************************************************************************************
                            qjobuvwidget.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <sstream>

#include "qjobuvwidget.h"
#include "qconvwidget.h"
#include "qirwidget.h"
#include "qjcdrawing.h"
#include "qplotspectrum.h"

#include "world.h"
#include "glvisor.h"
#include "molecule.h"
#include "frame.h"
#include "kryovisoroptical.h"

#include <QDockWidget>

QJobUVWidget::QJobUVWidget(const QString& file, QWidget* parent ) : QJobWidget (parent), m_file (file)
{
    m_world = new kryomol::World(this,kryomol::World::opticalvisor);
}

QJobUVWidget::~QJobUVWidget()
{

}

void QJobUVWidget::InitWidgets()
{
    this->World()->Initialize();
    this->setCentralWidget(m_world->Visor());
    InitCommonWidgets();



    QDockWidget* uvdock = new QDockWidget(this);
    uvdock->setAllowedAreas(Qt::RightDockWidgetArea);

    m_uvwidget = new QUVWidget (this->World(),m_file, m_tabwidget);
    m_tabwidget->addTab(m_uvwidget,"UV/VIS");

    uvdock->setWidget(m_tabwidget);


    this->addDockWidget(Qt::RightDockWidgetArea,uvdock);

    /*uvdock->show();
    for(auto it=m_dockwidgets.begin();it!=m_dockwidgets.end();++it)
    {
        this->tabifyDockWidget(uvdock,*it);
    }
    m_dockwidgets.insert(0,uvdock);*/
    //uvdock->raise();


    connect ( m_uvwidget,SIGNAL( Type ( QPlotSpectrum::SpectrumType ) ),this,SLOT( OnUVTypeChanged ( QPlotSpectrum::SpectrumType ) ) );
    connect ( m_uvwidget,SIGNAL( showspectrum ( bool ) ),this,SLOT ( OnShowUVSpectrum ( bool ) ) );
    connect ( m_uvwidget,SIGNAL( showelectricdipole(bool)),World()->Visor(),SLOT(OnShowElectricDipole(bool )));
    connect ( m_uvwidget,SIGNAL( showmagneticdipole(bool)),World()->Visor(),SLOT(OnShowMagneticDipole(bool )));
    connect ( m_uvwidget,SIGNAL( showvelocitydipole(bool)),World()->Visor(),SLOT(OnShowVelocityDipole(bool )));
    connect ( m_uvwidget,SIGNAL( setactivetransition(int)),World()->Visor(),SLOT(OnSetActiveTransition(int )));

    //this->World()->Visor()->Initialize();
    //Get all the lines
    std::vector< std::vector<Spectralline> > linesets;
    std::vector< std::vector< std::vector< kryomol::TransitionChange> > > tcsets;
    for( const auto& f : World()->CurrentMolecule()->Frames() )
    {
        linesets.push_back(f.GetSpectralLines());
        tcsets.push_back(f.TransitionChanges());
    }
    m_uvwidget->SetLines(linesets);
    m_uvwidget->SetTransitionChanges(tcsets);
    m_uvwidget->InitTable(World()->CurrentMolecule()->CurrentFrameIndex());

    connect(World(),SIGNAL(currentFrame(size_t)),this,SLOT(OnFrameChanged(size_t)));



}

void QJobUVWidget::OnShowUVSpectrum ( bool bshow )
{

    if ( m_irwidget )
    {
        bshow ? m_irwidget->show() : m_irwidget->hide();
    }
    if ( !bshow ) return;

    QDockWidget* irdock = new QDockWidget(this);
    m_irwidget = new QIRWidget(this->World(),this);
    irdock->setWidget(m_irwidget);
    irdock->setAllowedAreas(Qt::BottomDockWidgetArea);

    QPlotSpectrum* jc= m_irwidget->GetSpectrum();
    //jc->SetType(QPlotSpectrum::UV);
    connect ( m_uvwidget,SIGNAL ( data ( const std::vector<fidarray>*,const fidarray*,float,float,float,QPlotSpectrum::SpectrumType) ),jc,SLOT ( SetData ( const std::vector<fidarray>*,const fidarray*,float,float,float,QPlotSpectrum::SpectrumType) ) );
    std::vector<QColor> colors;
    for( const auto& f : this->World()->CurrentMolecule()->Frames() )
    {
        float h,s,l;
        f.GetColor(h,s,l);
        colors.push_back(QColor::fromHslF(h,s,l));
    }
    jc->SetColors(colors);
    jc->SetData ( m_uvwidget->GetData(),m_uvwidget->GetTotalData(),m_uvwidget->Max(),m_uvwidget->Min(), m_uvwidget->Shift(),QPlotSpectrum::UV);
    //set the colors


}


void QJobUVWidget::OnUVTypeChanged ( QPlotSpectrum::SpectrumType t )
{

#warning should not be QUVWidget?
    QList<QIRWidget*> childlist= this->findChildren<QIRWidget*>();
    if ( childlist.empty() ) return;
    QIRWidget* w=childlist.first();
    if ( ! w ) return;

    w->GetSpectrum()->SetType(t);
}

void QJobUVWidget::OnFrameChanged(size_t findex)
{
    m_uvwidget->InitTable(findex);
}
