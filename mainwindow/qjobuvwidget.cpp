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

QJobUVWidget::QJobUVWidget(const QString& file, kryomol::World* world, QWidget* parent ) : QJobWidget (world, parent), m_file (file)
{
   Init();
}

QJobUVWidget::~QJobUVWidget()
{
}

void QJobUVWidget::Init()
{
    kryomol::World* world = GetWorld();

    m_uvwidget = new QUVWidget (world,m_file, this);

    connect ( m_uvwidget,SIGNAL( Type ( QPlotSpectrum::SpectrumType ) ),this,SLOT( OnUVTypeChanged ( QPlotSpectrum::SpectrumType ) ) );
    connect ( m_uvwidget,SIGNAL( showspectrum ( bool ) ),this,SLOT ( OnShowUVSpectrum ( bool ) ) );
    connect ( m_uvwidget,SIGNAL( showelectricdipole(bool)),world->Visor(),SLOT(OnShowElectricDipole(bool )));
    connect ( m_uvwidget,SIGNAL( showmagneticdipole(bool)),world->Visor(),SLOT(OnShowMagneticDipole(bool )));
    connect ( m_uvwidget,SIGNAL( showvelocitydipole(bool)),world->Visor(),SLOT(OnShowVelocityDipole(bool )));
    connect ( m_uvwidget,SIGNAL( setactivetransition(int)),world->Visor(),SLOT(OnSetActiveTransition(int )));

    world->Visor()->Initialize();
    //Get all the lines
    std::vector< std::vector<Spectralline> > linesets;
    std::vector< std::vector< std::vector< kryomol::TransitionChange> > > tcsets;
    for( const auto& f : world->CurrentMolecule()->Frames() )
    {
        linesets.push_back(f.GetSpectralLines());
        tcsets.push_back(f.TransitionChanges());
    }
    m_uvwidget->SetLines(linesets);
    m_uvwidget->SetTransitionChanges(tcsets);
    m_uvwidget->InitTable(world->CurrentMolecule()->CurrentFrameIndex());

    //Add the visor and the FreqWidget to the splitter
    //this->addWidget(world->Visor());
    this->addWidget(m_uvwidget);

    SetWorld(world);

    connect(world,SIGNAL(currentFrame(size_t)),this,SLOT(OnFrameChanged(size_t)));

}

void QJobUVWidget::OnShowUVSpectrum ( bool bshow )
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

    QList<QUVWidget*> chw= this->findChildren<QUVWidget*>();
    QUVWidget* fw= chw.first();


    QPlotSpectrum* jc= ir->GetSpectrum();
    //jc->SetType(QPlotSpectrum::UV);
    connect ( fw,SIGNAL ( data ( const std::vector<fidarray>*,const fidarray*,float,float,float,QPlotSpectrum::SpectrumType) ),jc,SLOT ( SetData ( const std::vector<fidarray>*,const fidarray*,float,float,float,QPlotSpectrum::SpectrumType) ) );
    std::vector<QColor> colors;
    for( const auto& f : this->GetWorld()->CurrentMolecule()->Frames() )
    {
        float h,s,l;
        f.GetColor(h,s,l);
        colors.push_back(QColor::fromHslF(h,s,l));
    }
    jc->SetColors(colors);
    jc->SetData ( fw->GetData(),fw->GetTotalData(),fw->Max(),fw->Min(), fw->Shift(),QPlotSpectrum::UV);
    //set the colors



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
