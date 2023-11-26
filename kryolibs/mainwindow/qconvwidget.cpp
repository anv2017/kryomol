/*****************************************************************************************
                            qconvwidget.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <iostream>

#include <qslider.h>

#include <QLabel>
#include <QLayout>
#include <QPushButton>
#include <QFrame>
#include <QCheckBox>
#include <QGroupBox>
#include <qwt_plot.h>

#include "qconvwidget.h"
//#include "qryoplot.h"
#include "quantumplot.h"
#include "qdoubleeditbox.h"


Convergence::Convergence ( size_t size )
{
    m_energies= new double[size];
    m_S2= new double[size];
    m_rmsforces= new double[size];
    m_maximumforces= new double[size];
    m_rmsdisplacements= new double[size];
    m_maximumdisplacements= new double[size];
    m_size=size;
}

Convergence::~Convergence()
{
    delete [] m_energies;
    delete [] m_S2;
    delete [] m_rmsforces;
    delete [] m_maximumforces;
    delete [] m_rmsdisplacements;
    delete [] m_maximumdisplacements;

}

double Convergence::MaxEnergy()
{
    return GetMaxArray ( m_energies );

}

double Convergence::MinEnergy()
{
    return GetMinArray ( m_energies );
}

double Convergence::GetMaxArray ( double* array )
{
    double max=-1e25;
    for ( size_t i=0;i< m_size;i++ )
    {
        if ( array[i] > max )
            max=array[i];
    }

    return max;
}

double Convergence::GetMinArray ( double* array )
{
    double min=1e25;
    for ( size_t i=0;i< m_size;i++ )
    {
        if ( array[i] < min )
            min=array[i];
    }
    return min;
}

std::ostream& operator << ( std::ostream& s, Convergence& conv )
{
    s << "Energy" << "\t" << " RMS Force" << "\t" << "Max.Force" << "\t" << "RMS Disp." << "\t" << "MaxDisp" << std::endl;

    double* energies=conv.Energies();
    double* rmsforces=conv.RMSForces();
    double* maxforces=conv.MaximumForces();
    double* rmsdisp=conv.RMSDisplacements();
    double* maxdisp=conv.MaximumDisplacements();

    for ( size_t i=0;i < conv.size();i++ )
    {
        s << energies[i] << "\t" << maxforces[i] << "\t" << rmsforces[i] <<"\t"
          << maxdisp[i] <<"\t" << rmsdisp[i] << std::endl;

    }
    s << std::endl;
    s << energies[0] << std::endl;
    return s;
}

QConvWidget::QConvWidget (QWidget *parent,bool gradients /*=true*/, bool displacements /*=true*/)
    : QWidget ( parent ) , m_gradient(NULL), m_displacement(NULL)
{
    m_bshowforces = false;

    QVBoxLayout* lay = new QVBoxLayout ( this );

    m_energy = new QuantumPlot ( this );
    lay->addWidget ( m_energy );
    lay->insertSpacing(1,50);

    if ( gradients )
    {
        m_gradient = new QuantumPlot ( this );

        lay->addWidget ( m_gradient );
        lay->insertSpacing(3,50);
    }

    if ( displacements)
    {
        m_displacement = new QuantumPlot ( this );

        lay->addWidget ( m_displacement  );
        lay->insertSpacing(5,50);
    }

    m_limitsbox = new QDoubleEditBox ( 0,0,this );
    connect ( m_limitsbox,SIGNAL ( limits ( int,int ) ),this,SLOT ( OnChangeLimits ( int,int ) ) );
    m_limitsbox->setWhatsThis(tr(""));


    lay->addWidget ( m_limitsbox, 0, 0 );

    if ( gradients )
    {
        QGroupBox* box = new QGroupBox(this);
        box->setTitle("Forces");
        box->setWhatsThis(tr(""));

        QGridLayout* hlay = new QGridLayout(this);

        m_showforce = new QCheckBox ("Show",this);
        hlay->addWidget(m_showforce,0,0);
        connect (m_showforce, SIGNAL(clicked()),this,SLOT (OnShowForces()));

        m_scaleslider = new QSlider (Qt::Horizontal,this );
        m_scaleslider->setTickPosition(QSlider::TicksBelow);
        hlay->addWidget ( m_scaleslider,0,2 );
        QLabel* slidertext= new QLabel ( this );
        slidertext->setText ( "Scale" );
        hlay->addWidget ( slidertext,1,2,Qt::AlignCenter );
        connect ( m_scaleslider,SIGNAL ( valueChanged ( int ) ),this,SLOT ( OnForceSliderChanged ( int ) ) );

        m_scaleslider->setValue(m_scaleslider->maximum()/2);

        hlay->setSizeConstraint(QLayout::SetMinimumSize);

        box->setLayout(hlay);

        lay->addWidget(box);

    }
    m_xaxis=NULL;
    m_convdata=NULL;

    connect ( m_energy,SIGNAL ( selectedPoint ( size_t ) ),this,SIGNAL ( selectedPoint ( size_t ) ) );
    connect ( m_energy,SIGNAL ( selectedPoint ( size_t ) ),m_gradient,SLOT ( OnSelectedPoint ( size_t ) ) );
    connect ( m_energy,SIGNAL ( selectedPoint ( size_t ) ),m_displacement,SLOT ( OnSelectedPoint ( size_t ) ) );

    if (gradients )
    {
        connect ( m_gradient,SIGNAL ( selectedPoint ( size_t ) ),this,SIGNAL ( selectedPoint ( size_t ) ) );
        connect ( m_gradient,SIGNAL ( selectedPoint ( size_t ) ),m_energy,SLOT ( OnSelectedPoint ( size_t ) ) );
        connect ( m_gradient,SIGNAL ( selectedPoint ( size_t ) ),m_displacement,SLOT ( OnSelectedPoint ( size_t ) ) );
    }

    if ( displacements )
    {
        connect ( m_displacement,SIGNAL ( selectedPoint ( size_t ) ),this,SIGNAL ( selectedPoint ( size_t ) ) );
        connect ( m_displacement,SIGNAL ( selectedPoint ( size_t ) ),m_energy,SLOT ( OnSelectedPoint ( size_t ) ) );
        connect ( m_displacement,SIGNAL ( selectedPoint ( size_t ) ),m_gradient,SLOT ( OnSelectedPoint ( size_t ) ) );
    }

}


QConvWidget::~QConvWidget()
{

    if ( m_xaxis ) delete [] m_xaxis;
    if ( m_convdata ) delete m_convdata;
}

void QConvWidget::SetNData ( size_t size )
{

    m_xaxis=new double[size];
    m_ndata=size;

    for ( size_t  i=0;i< size;i++ )
    {
        m_xaxis[i]=i+1;
    }

    m_convdata= new Convergence ( m_ndata );
    m_limitsbox->SetLimits ( 1,m_ndata );


}

void QConvWidget::SetupCurves()
{
    //Take into account incomplete minimizations
    double* forces=m_convdata->RMSForces();
    int data=m_ndata;
    if ( forces[m_ndata-1] < 0 )
    {
        data=m_ndata-1;
    }
    EnergyCurve* envc= new EnergyCurve ( m_energy,"energies" );
    envc->setExtraLabel("E");
    envc->setRawSamples( m_xaxis,m_convdata->Energies(),data  );
    envc->setExtraData ( m_convdata->S2() );
    envc->setPen ( QPen ( QColor ( 255,0,0 ) ) );
    envc->attach ( m_energy );
    m_energy->setAxisTitle ( QwtPlot::yLeft,"Energy" );
    m_energy->setupStepMarkers();
    m_energy->replot();

    if ( m_gradient )
    {
        QuantumCurve* mxfcv= new QuantumCurve ( m_gradient,"maxforces" );
        mxfcv->setExtraLabel("Max. Force");
        QuantumCurve* rmsfcv= new QuantumCurve ( m_gradient,"rmsforces" );
        rmsfcv->setExtraLabel("RMS Force");

        rmsfcv->setRawSamples ( m_xaxis,m_convdata->RMSForces(),data );
        mxfcv->setRawSamples ( m_xaxis,m_convdata->MaximumForces(),data );
        rmsfcv->setPen ( QPen ( QColor ( 255,0,0 ) ) );
        mxfcv->setPen ( QPen ( QColor ( 0,0,255 ) ) );
        rmsfcv->attach ( m_gradient );
        mxfcv->attach ( m_gradient );
        m_gradient->setAxisTitle ( QwtPlot::yLeft,"Forces" );
        QwtPlotMarker* maxfmarker= new QwtPlotMarker();
        maxfmarker->setLineStyle ( QwtPlotMarker::HLine );
        maxfmarker->setYValue ( m_threshold.maxforce );
        maxfmarker->setLinePen ( QPen ( Qt::green,0,Qt::DotLine ) );
        QwtPlotMarker* rmsfmarker=new QwtPlotMarker();
        rmsfmarker->setLineStyle ( QwtPlotMarker::HLine );
        rmsfmarker->setYValue ( m_threshold.rmsforce );
        rmsfmarker->setLinePen ( QPen ( Qt::yellow,0,Qt::DotLine ) );
        maxfmarker->attach ( m_gradient );
        rmsfmarker->attach ( m_gradient );
        m_gradient->setupStepMarkers();
        m_gradient->replot();
    }
    if ( m_displacement )

    {
        QuantumCurve* rmsdisp = new QuantumCurve ( m_displacement,"rmsdisplacements" );
        rmsdisp->setExtraLabel("RMS Disp.");
        QuantumCurve* maxdisp = new QuantumCurve ( m_displacement,"maxdisplacements" );
        maxdisp->setExtraLabel("Max. Disp.");

        rmsdisp->setRawSamples ( m_xaxis,m_convdata->RMSDisplacements(),data );
        maxdisp->setRawSamples ( m_xaxis,m_convdata->MaximumDisplacements(),data );
        rmsdisp->setPen ( QPen ( QColor ( 255,0,0 ) ) );
        maxdisp->setPen ( QPen ( QColor ( 0,0,255 ) ) );

        rmsdisp->attach ( m_displacement );
        maxdisp->attach ( m_displacement );

        m_displacement->setAxisTitle ( QwtPlot::yLeft,"Displacements" );
        QwtPlotMarker* maxdmarker=new QwtPlotMarker();;
        maxdmarker->setLineStyle ( QwtPlotMarker::HLine );
        maxdmarker->setYValue ( m_threshold.maxdisplacement );
        maxdmarker->setLinePen ( QPen ( Qt::green,0,Qt::DotLine ) );
        QwtPlotMarker* rmsdmarker=new QwtPlotMarker();
        rmsdmarker->setLineStyle ( QwtPlotMarker::HLine );
        rmsdmarker->setYValue ( m_threshold.rmsdisplacement );
        rmsdmarker->setLinePen ( QPen ( Qt::yellow,0,Qt::DotLine ) );
        maxdmarker->attach ( m_displacement );
        rmsdmarker->attach ( m_displacement );
        m_displacement->setupStepMarkers();
        m_displacement->replot();
    }

}


void QConvWidget::OnChangeLimits ( int first, int last )
{

    double max=-1e30;
    double min=1e30;
    for ( int i=first-1;i<last;++i )
    {
        if ( m_convdata->Energies() [i] > max ) max= ( m_convdata->Energies() [i] );
        if ( m_convdata->Energies() [i] < min ) min= ( m_convdata->Energies() [i] );
    }
    m_energy->ChangeLimits ( first,last,min,max );
    if ( m_gradient )
    {
        max=-1e30;
        min=1e30;
        for ( int i=first-1;i<last;++i )
        {
            if ( m_convdata->MaximumForces() [i] > max ) max= ( m_convdata->MaximumForces() [i] );
            if ( m_convdata->RMSForces() [i] > max ) max= ( m_convdata->RMSForces() [i] );
            if ( m_convdata->MaximumForces() [i] < min ) min= ( m_convdata->MaximumForces() [i] );
            if ( m_convdata->RMSForces() [i] < min ) min= ( m_convdata->RMSForces() [i] );

        }
        m_gradient->ChangeLimits(first,last,min,max);
    }
    if ( m_displacement )
    {
        max=-1e30;
        min=1e30;
        for ( int i=first-1;i<last;++i )
        {
            if ( m_convdata->MaximumDisplacements() [i] > max ) max= ( m_convdata->MaximumDisplacements() [i] );
            if ( m_convdata->RMSDisplacements() [i] > max ) max= ( m_convdata->RMSDisplacements() [i] );
            if ( m_convdata->MaximumDisplacements() [i] < min ) min= ( m_convdata->MaximumDisplacements() [i] );
            if ( m_convdata->RMSDisplacements() [i] < min ) min= ( m_convdata->RMSDisplacements() [i] );

        }

        m_displacement->ChangeLimits(first,last,min,max);
    }
}


void QConvWidget::OnSelectedPoint ( size_t point )
{
    m_energy->OnSelectedPoint ( point );
    if ( m_gradient )
        m_gradient->OnSelectedPoint ( point );
    if ( m_displacement )
        m_displacement->OnSelectedPoint ( point );
}


void QConvWidget::OnShowForces()
{
    m_bshowforces =! m_bshowforces;
    emit showforces(m_bshowforces);
}


void QConvWidget::OnForceSliderChanged ( int v )
{
    emit forcescale ( sqrt(v) );
}
