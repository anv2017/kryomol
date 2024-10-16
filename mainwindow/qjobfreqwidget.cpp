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
#include "qplotspectrum.h"
//#include "qjcdrawing.h"

#include "world.h"
#include "kryovisor.h"

QJobFreqWidget::QJobFreqWidget(const QString& file, kryomol::World* world, QWidget* parent ) : QJobWidget (world, parent), m_file (file)
{
    Init();
}

QJobFreqWidget::~QJobFreqWidget()
{
}


void QJobFreqWidget::Init()
{
    kryomol::World* world = GetWorld();

    //Create the FreqWidget
    m_freqwidget = new QFreqWidget (world,m_file, this);

    connect ( m_freqwidget,SIGNAL ( Type ( QPlotSpectrum::SpectrumType ) ),this,SLOT ( OnIRTypeChanged ( QPlotSpectrum::SpectrumType ) ) );

    //Initialize the visor and actions: spectrum,...
    connect ( m_freqwidget,SIGNAL ( mode ( int,int ) ),world->Visor(),SLOT ( SetMode ( int,int ) ) );
    connect ( m_freqwidget,SIGNAL ( distort ( int ) ),world->Visor(),SLOT ( OnDistortion ( int ) ) );
    connect ( m_freqwidget,SIGNAL ( reset() ),world->Visor(),SLOT ( OnStopAnimation() ) );
    connect ( m_freqwidget,SIGNAL ( showspectrum ( bool ) ),this,SLOT ( OnShowSpectrum ( bool ) ) );
    world->Visor()->Initialize();

    m_freqwidget->InitFrequencies();
    m_freqwidget->InitTable(0);

    //Add the visor and the FreqWidget to the splitter
    this->addWidget(world->Visor());
    this->addWidget(m_freqwidget);

    SetWorld(world);
}


void QJobFreqWidget::OnShowSpectrum ( bool bshow )
{
  //if already exist a irwidget object simply show it, otherwise create it
  if ( bshow )
  {
    QList<QIRWidget*> childlist= this->findChildren<QIRWidget*>(); ;
    if ( !childlist.empty() )
    {

      for ( QList<QIRWidget*>::iterator ch=childlist.begin();ch!=childlist.end();++ch)
         (*ch)->show();
      return;
    }


    //add a vertical splitter to the current splitter
    QSplitter* vsplit =new QSplitter ( Qt::Vertical, this );

    //now get a pointer to the first element of currennt split
    kryomol::KryoVisor* w= this->findChild<kryomol::KryoVisor*>();


    //and now reparent this widget which should be the GLVisor
    w->setParent(vsplit,0);
    w->move(QPoint(0,0));

    w->setSizePolicy ( QSizePolicy::Expanding,QSizePolicy::Expanding );
    vsplit->addWidget(w);
    //add the Spectrum drawing widget
    QIRWidget* ir= new QIRWidget ( this->GetWorld(),vsplit );
    vsplit->addWidget(ir);
    vsplit->setStretchFactor ( vsplit->indexOf(ir), 0);

    QList<QFreqWidget*> chw= this->findChildren<QFreqWidget*>();
    QFreqWidget* fw= chw.first();


    QPlotSpectrum* jc= ir->GetSpectrum();
    connect ( fw,SIGNAL ( data ( fidarray*,float,float, float ) ),jc,SLOT ( SetData ( fidarray*,float,float,float ) ) );
    jc->SetData ( &fw->GetData(),&fw->GetTotalData(),fw->Max(),fw->Min(), fw->Shift(),QPlotSpectrum::IR);

    //Change the size of the widgets in the splitter
    QList<int> list = vsplit->sizes();
    list[0]=fw->size().height()/2;
    list[1]=fw->size().height()/2;
    vsplit->setSizes(list);

    //and finally move the vertical splitter to the first position
    this->insertWidget(0, vsplit);

    list = this->sizes();
    list[0]=this->size().width()-this->size().width()/4;
    list[1]=this->size().width()/4;
    this->setSizes(list);

  }
  else
  {
    //find a QJCDrawing and hide it
    QList<QIRWidget*> childlist=this->findChildren<QIRWidget*>();

    for ( QList<QIRWidget*>::iterator ch=childlist.begin();ch!=childlist.end();++ch)
      (*ch)->hide();
  }
}


void QJobFreqWidget::OnIRTypeChanged ( QPlotSpectrum::SpectrumType t )
{

#warning should not be QUVWidget?
  QList<QIRWidget*> childlist= this->findChildren<QIRWidget*>();
  if ( childlist.empty() ) return;
  QIRWidget* w=childlist.first();
  if ( ! w ) return;

  w->GetSpectrum()->SetType(t);
}

