/*****************************************************************************************
                            qirwidget.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <qlayout.h>
#include <qframe.h>
#include <qpushbutton.h>

#include "qirwidget.h"
#include "qjcdrawing.h"
#include "qplotspectrum.h"
//Added by qt3to4:
#include <QVBoxLayout>
#include <QPixmap>

#include "confmanager.h"


QIRWidget::QIRWidget(kryomol::World* w, QWidget* parent ) :  QWidget(parent)
{  
    setupUi(this);
    QHBoxLayout* hlay= new QHBoxLayout(_spectrumFrame);

    m_spectrum= new QPlotSpectrum(_spectrumFrame);
    m_spectrum->setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Expanding);
    m_confmanager = new ConfManager(w,this);
    m_confmanager->setVisible(false);
    hlay->addWidget(m_spectrum);
    hlay->addWidget(m_confmanager);



    InitButtons();

}

void QIRWidget::InitButtons()
{
    connect(_printButton,SIGNAL(clicked()),m_spectrum,SLOT(OnPrint()));
    connect(_plusButton,SIGNAL(clicked()),m_spectrum,SLOT(OnIncrease()));
    connect(_minusButton,SIGNAL(clicked()),m_spectrum,SLOT(OnDecrease()));
    connect(_zoomButton,SIGNAL(toggled(bool)),m_spectrum,SLOT(OnZoom(bool)));
    connect(_fullButton,SIGNAL(clicked()),m_spectrum,SLOT(ResetZoom()));
    connect(_confButton,SIGNAL(toggled(bool)),this,SLOT(OnShowConformers(bool)));

    QAction* bitmapAction = new QAction( tr("Bitmap"),this);
    connect ( bitmapAction,SIGNAL ( activated() ),m_spectrum, SLOT ( OnBitmapPicture() ) );
    _pictureButton->addAction(bitmapAction);
    QAction* SVGAction = new QAction( tr("SVG"),this);
    connect ( SVGAction,SIGNAL ( activated() ),m_spectrum, SLOT ( OnSVGPicture() ) );
    _pictureButton->addAction(SVGAction);
    _pictureButton->setPopupMode(QToolButton::InstantPopup);
}

QIRWidget::~QIRWidget()
{}

void QIRWidget::OnShowConformers(bool b)
{
    m_confmanager->setVisible(b);
}


