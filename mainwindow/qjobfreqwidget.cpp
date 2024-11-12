/*****************************************************************************************
                            qjobfreqwidget.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "qjobfreqwidget.h"
#include "qjobwidget.h"
#include "qirwidget.h"
#include "qfreqwidget.h"
#include "qplotspectrum.h"
//#include "qjcdrawing.h"

#include "world.h"
#include "kryovisor.h"

#include <QDockWidget>

QJobFreqWidget::QJobFreqWidget(const QString& file, QWidget* parent ) : QJobWidget (parent),
    m_file (file), m_irwidget(nullptr)
{
       m_world = new kryomol::World(this,kryomol::World::freqvisor);
}

QJobFreqWidget::~QJobFreqWidget()
{
}


void QJobFreqWidget::InitWidgets()
{
    this->World()->Initialize();
    this->setCentralWidget(m_world->Visor());
    InitCommonWidgets();


    //Create the FreqWidget
    QDockWidget* fdock = new QDockWidget(this);
    fdock->setAllowedAreas(Qt::RightDockWidgetArea);
    m_freqwidget = new QFreqWidget (this->World(),m_file, m_tabwidget);
    m_tabwidget->addTab(m_freqwidget,"Freq");
    fdock->setWidget(m_tabwidget);
    this->addDockWidget(Qt::RightDockWidgetArea,fdock);

    connect ( m_freqwidget,SIGNAL ( Type ( QPlotSpectrum::SpectrumType ) ),this,SLOT ( OnIRTypeChanged ( QPlotSpectrum::SpectrumType ) ) );

    //Initialize the visor and actions: spectrum,...
    connect ( m_freqwidget,SIGNAL ( mode ( int,int ) ),World()->Visor(),SLOT ( SetMode ( int,int ) ) );
    connect ( m_freqwidget,SIGNAL ( distort ( int ) ),World()->Visor(),SLOT ( OnDistortion ( int ) ) );
    connect ( m_freqwidget,SIGNAL ( reset() ),World()->Visor(),SLOT ( OnStopAnimation() ) );
    connect ( m_freqwidget,SIGNAL ( showspectrum ( bool ) ),this,SLOT ( OnShowSpectrum ( bool ) ) );
    World()->Visor()->Initialize();

    m_freqwidget->InitFrequencies();
    m_freqwidget->InitTable(0);

}


void QJobFreqWidget::OnShowSpectrum ( bool bshow )
{

    if ( m_irwidget )
    {
        bshow ?  m_irwidget->show() : m_irwidget->hide();
        return;
    }

    if ( !bshow ) return;

    //Create the IRWidget
    QDockWidget* irdock = new QDockWidget(this);
    this->addDockWidget(Qt::BottomDockWidgetArea,irdock);
    m_irwidget= new QIRWidget ( this->World(),irdock );
    irdock->setWidget(m_irwidget);
    irdock->setAllowedAreas(Qt::BottomDockWidgetArea);

    QPlotSpectrum* jc= m_irwidget->GetSpectrum();
    connect ( m_freqwidget,SIGNAL ( data ( fidarray*,float,float, float ) ),jc,SLOT ( SetData ( fidarray*,float,float,float ) ) );
    std::vector<QColor> colors;
    for( const auto& f : this->World()->CurrentMolecule()->Frames() )
    {
        float h,s,l;
        f.GetColor(h,s,l);
        colors.push_back(QColor::fromHslF(h,s,l));
    }
    jc->SetColors(colors);
    jc->SetData ( &m_freqwidget->GetData(),&m_freqwidget->GetTotalData(),m_freqwidget->Max(),m_freqwidget->Min(), m_freqwidget->Shift(),QPlotSpectrum::IR);


}


void QJobFreqWidget::OnIRTypeChanged ( QPlotSpectrum::SpectrumType t )
{
  m_irwidget->GetSpectrum()->SetType(t);
}

