/*****************************************************************************************
                            qorbitalwidget.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/


#include "qorbitalwidget.h"
#include "ui_qorbitalwidgetbase.h"
#include "qnumlistviewitem.h"
#include "density.h"
#include "qplotspectrogram.h"
#include "renderorbitals.h"
#include "molecule.h"
#include "frame.h"

#include "QButtonGroup"
#include <QTime>

#include <sstream>
#include <iostream>


QOrbitalWidget::QOrbitalWidget(QWidget *parent) :
    QWidget(parent)
{
    setupUi(this);

    Init();
}

QOrbitalWidget::~QOrbitalWidget()
{

}

void QOrbitalWidget::Init()
{
    m_world=nullptr;
    m_bshowhomo = false;
    m_bshowlumo = false;
    m_bshowtotaldensity = false;
    m_bshowspindensity = false;
    m_bonshow = false;
    m_bshowcontours = false;
    m_orbital = 0;
    m_isovalue = 0.0032;
    m_gridresolution = 0.2;
    m_threshold = 0.001;

    _sliderIsovalue->setOrientation(Qt::Horizontal);
    _sliderIsovalue->setTickPosition(QSlider::TicksBelow);
    _sliderIsovalue->setTickInterval(1);
    _sliderIsovalue->setMaximum(10);
    _sliderIsovalue->setMinimum(0);
    _sliderIsovalue->setDoubleMaximum(1.0);
    _sliderIsovalue->setDoubleMinimum(0.00001);
    _sliderIsovalue->setValue(5);
    _sliderIsovalue->setTracking(false);
    connect(_sliderIsovalue,SIGNAL(valueChanged(int)),_sliderIsovalue,SLOT(OnValueChanged(int)));
    connect(_sliderIsovalue,SIGNAL(valueChanged(double)),this,SLOT(OnIsovalueSliderChange(double)));

    connect(_sliderTransparence,SIGNAL(valueChanged(int)),this,SLOT(OnTransparenceSliderChange(int)));

    connect(_buttonTotalDensity,SIGNAL (clicked()),this,SLOT(OnShowTotalDensity()));

    connect(_buttonShowContours,SIGNAL(clicked()),this,SLOT(OnShowContoursFrame()));

    _labelIsovalue->show();
    _sliderIsovalue->show();
    _labelTransparence->show();
    _sliderTransparence->show();
    _labelGridResolution->hide();
    _sliderGridResolution->hide();
    _labelThreshold->hide();
    _thresholdSpinBox->hide();
    _labelOrbitals->hide();
    _labelOrbitals2->hide();
    _listOrbitals->hide();
    _listOrbitals2->hide();
    _buttonHomo->hide();
    _buttonLumo->hide();
    _buttonTotalDensity->show();
    _buttonSpinDensity->hide();
    _buttonShowContours->show();
    _contoursGroupBox->hide();


    QVBoxLayout* vlay= new QVBoxLayout(_frameSpectrogram);

    m_plotspectrogram= new QPlotSpectrogram(&m_render.DensityData(),_frameSpectrogram);
    m_plotspectrogram->setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Expanding);

    vlay->addWidget(m_plotspectrogram);

    _buttonSpectrogram->setCheckable(true);
    _buttonContour->setCheckable(true);

    QButtonGroup* _buttonGroupAxis = new QButtonGroup(_contoursGroupBox);
    _buttonGroupAxis->addButton(_checkBoxXY);
    _buttonGroupAxis->addButton(_checkBoxXZ);
    _buttonGroupAxis->addButton(_checkBoxYZ);
    _buttonGroupAxis->setExclusive(true);
    _checkBoxXY->setCheckable(true);
    _checkBoxXZ->setCheckable(true);
    _checkBoxYZ->setCheckable(true);
    _checkBoxXY->setChecked(true);
    _checkBoxXZ->setChecked(false);
    _checkBoxYZ->setChecked(false);
    _spinBoxXY->setEnabled(true);
    _spinBoxXZ->setEnabled(false);
    _spinBoxYZ->setEnabled(false);

#if QT_VERSION >= 0x040000
    _buttonSpectrogram->setChecked(true);
    _buttonContour->setChecked(false);
#else
    _buttonSpectrogram->setOn(true);
    _buttonContour->setOn(false);
#endif

    connect(_buttonSpectrogram, SIGNAL(toggled(bool)), m_plotspectrogram, SLOT(showSpectrogram(bool)));
    connect(_buttonContour, SIGNAL(toggled(bool)), m_plotspectrogram, SLOT(showContour(bool)));
    connect(_buttonPrint, SIGNAL(clicked()), m_plotspectrogram, SLOT(printSpectrogram()) );

    connect(_checkBoxXY, SIGNAL(toggled(bool)), this, SLOT(OnChangeAxisPlane(bool)));
    connect(_checkBoxXZ, SIGNAL(toggled(bool)), this, SLOT(OnChangeAxisPlane(bool)));
    connect(_checkBoxYZ, SIGNAL(toggled(bool)), this, SLOT(OnChangeAxisPlane(bool)));
    connect(_spinBoxXY, SIGNAL(valueChanged(double)), this, SLOT(OnChangeAxisValue(double)));
    connect(_spinBoxXZ, SIGNAL(valueChanged(double)), this, SLOT(OnChangeAxisValue(double)));
    connect(_spinBoxYZ, SIGNAL(valueChanged(double)), this, SLOT(OnChangeAxisValue(double)));
}

void QOrbitalWidget::SetRenderOrbitals(kryomol::RenderOrbitals render)
{
    _labelIsovalue->show();
    _sliderIsovalue->show();
    _labelTransparence->show();
    _sliderTransparence->show();
    _labelGridResolution->hide();
    _sliderGridResolution->hide();
    _labelThreshold->hide();
    _thresholdSpinBox->hide();
    _labelOrbitals->hide();
    _labelOrbitals2->hide();
    _listOrbitals->hide();
    _listOrbitals2->hide();
    _buttonHomo->hide();
    _buttonLumo->hide();
    _buttonTotalDensity->show();
    _buttonSpinDensity->hide();
    _buttonShowContours->show();

    OffButtons(true);

    m_render=render;

    _spinBoxXY->setMaximum((m_render.DensityData().Nz()-1)*m_render.DensityData().Dz()/2);
    _spinBoxXY->setMinimum(-(m_render.DensityData().Nz()-1)*m_render.DensityData().Dz()/2);
    _spinBoxXY->setSingleStep(m_render.DensityData().Dz());
    _spinBoxXZ->setMaximum((m_render.DensityData().Ny()-1)*m_render.DensityData().Dy()/2);
    _spinBoxXZ->setMinimum(-(m_render.DensityData().Ny()-1)*m_render.DensityData().Dy()/2);
    _spinBoxXZ->setSingleStep(m_render.DensityData().Dy());
    _spinBoxYZ->setMaximum((m_render.DensityData().Nx()-1)*m_render.DensityData().Dx()/2);
    _spinBoxYZ->setMinimum(-(m_render.DensityData().Nx()-1)*m_render.DensityData().Dx()/2);
    _spinBoxYZ->setSingleStep(m_render.DensityData().Dx());

}

void QOrbitalWidget::ListOrbitals()
{
    disconnect(_listOrbitals, SIGNAL(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)),this,SLOT(OnOrbitalChange(QTreeWidgetItem*)));
    _listOrbitals->clear();
    connect(_listOrbitals, SIGNAL(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)),this,SLOT(OnOrbitalChange(QTreeWidgetItem*)));

    _sliderGridResolution->setTracking(false);
    _thresholdSpinBox->setMouseTracking(false);
    _thresholdSpinBox->setKeyboardTracking(false);
    _labelGridResolution->show();
    _sliderGridResolution->show();
    _labelThreshold->show();
    _thresholdSpinBox->show();

    m_render.SetBeta(m_beta);

    if (m_beta)
    {
        disconnect(_listOrbitals2, SIGNAL(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)),this,SLOT(OnBetaOrbitalChange(QTreeWidgetItem*)));
        _listOrbitals2->clear();
        connect(_listOrbitals2, SIGNAL(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)),this,SLOT(OnBetaOrbitalChange(QTreeWidgetItem*)));

        _labelOrbitals->setText("Alpha Molecular Orbitals");

        QStringList label;
        label << "MO" << "E(a.u.)" << "Occ.";
        _listOrbitals2->setColumnCount(3);
        _listOrbitals2->setHeaderLabels(label);

        std::vector<float> eigenvalues = m_render.OrbitalsData().BetaEigenvalues();

        for ( size_t i=0; i<eigenvalues.size(); ++i)
        {
            QString item, eigenvalue, occupation;

            item.sprintf ( "%d",static_cast<int> ( i+1 ) );
            eigenvalue.sprintf("%.5f",eigenvalues.at(i));
            occupation.sprintf("%.1f",0.0);

            QStringList sl;
            sl << item << eigenvalue << occupation;

            kryomol::QNumTreeWidgetItem* mitem= new kryomol::QNumTreeWidgetItem ( _listOrbitals2, sl );
            _listOrbitals2->addTopLevelItem(mitem);

        }

        connect(_listOrbitals2, SIGNAL(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)),this,SLOT(OnBetaOrbitalChange(QTreeWidgetItem*)));

        _listOrbitals2->setSortingEnabled(true);
        _listOrbitals2->sortItems(0,Qt::AscendingOrder);
        for( int i=0;i<3;++i)
            _listOrbitals2->resizeColumnToContents(i);

        _listOrbitals2->setEnabled(true);
        _labelOrbitals2->show();
        _listOrbitals2->show();

        _buttonSpinDensity->setEnabled(true);
        _buttonSpinDensity->show();
        connect(_buttonSpinDensity, SIGNAL (clicked()), this, SLOT (OnShowSpinDensity()));
    }

    QStringList label;
    label << "MO" << "E(a.u.)" << "Occ.";
    _listOrbitals->setColumnCount(3);
    _listOrbitals->setHeaderLabels(label);

    const std::vector<float>& eigenvalues = m_render.OrbitalsData().Eigenvalues();
    const std::vector<float>& occupations=m_render.OrbitalsData().Occupations();

    for ( size_t i=0; i<eigenvalues.size(); ++i)
    {
        QString item, eigenvalue, occupation;

        item.sprintf ( "%d",static_cast<int> ( i+1 ) );
        eigenvalue.sprintf("%.5f",eigenvalues.at(i));
        if ( !occupations.empty() )
            occupation.sprintf("%.1f",occupations.at(i));
        else occupation.sprintf("***");

        QStringList sl;
        sl << item << eigenvalue << occupation;

        kryomol::QNumTreeWidgetItem* mitem= new kryomol::QNumTreeWidgetItem ( _listOrbitals, sl );
        _listOrbitals->addTopLevelItem(mitem);

    }

    _listOrbitals->setSortingEnabled(true);
    _listOrbitals->sortItems(0,Qt::AscendingOrder);
    for( int i=0;i<3;++i)
        _listOrbitals->resizeColumnToContents(i);

    _buttonHomo->setEnabled(true);
    _buttonLumo->setEnabled(true);
    _listOrbitals->setEnabled(true);
    _thresholdSpinBox->setEnabled(true);

    _buttonHomo->show();
    _buttonLumo->show();
    _labelOrbitals->show();
    _listOrbitals->show();
    _thresholdSpinBox->show();

    connect(_listOrbitals, SIGNAL(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)),this,SLOT(OnOrbitalChange(QTreeWidgetItem*)));
    connect(_buttonHomo, SIGNAL(clicked()),this, SLOT(OnShowHomo()));
    connect(_buttonLumo, SIGNAL(clicked()),this, SLOT(OnShowLumo()));
    connect(_thresholdSpinBox, SIGNAL(valueChanged(double)),this, SLOT(OnThresholdAOSelectorChange(double)));
    connect(_sliderGridResolution, SIGNAL(valueChanged(int)),this, SLOT(OnGridResolutionChange(int)));
}

void QOrbitalWidget::OnIsovalueSliderChange(double isovalue)
{
    m_isovalue = isovalue;
    m_render.DensityData().SetIsovalue(isovalue);

    if (m_render.DensityData().ExistsDensityData())
        m_render.DensityData().RenderDensityData();

    if (m_bonshow)
        emit(drawDensity(true));
}


void QOrbitalWidget::OnTransparenceSliderChange(int value)
{
    float transparence = (float) value/10;
    emit(transparenceChange(transparence));
}

void QOrbitalWidget::OnThresholdAOSelectorChange(double v)
{
    m_threshold = v;
    m_render.SetThresholdAOSelector(v);
    m_render.CleanMatrices(m_bshowhomo, m_bshowlumo, m_bshowtotaldensity, m_bshowspindensity, m_orbital, m_betaorbital);

    emit drawDensity(m_bonshow);

    if ((m_bonshow)&&(!_contoursGroupBox->isHidden()))
        m_plotspectrogram->FillSpectrogram();
}


void QOrbitalWidget::OnGridResolutionChange(int v)
{
    float resolution;
    if (v<=4)
        resolution = 1.0-0.2*v;
    if (v>4)
        resolution = 0.1-0.02*(v-5);

    m_gridresolution = resolution;
    m_render.SetGridResolution(resolution);
    m_render.CleanMatrices(m_bshowhomo, m_bshowlumo, m_bshowtotaldensity, m_bshowspindensity, m_orbital, m_betaorbital);

    emit drawDensity(m_bonshow);

    if ((m_bonshow)&&(!_contoursGroupBox->isHidden()))
        m_plotspectrogram->FillSpectrogram();

}



void QOrbitalWidget::OnShowHomo()
{
    m_bshowhomo = !m_bshowhomo;
    m_bshowtotaldensity = false;
    m_bshowspindensity = false;
    m_bshowlumo = false;
    m_orbital = 0;//m_render.OrbitalsData().Homo();
    m_betaorbital = 0;

    emit offtransitions(true);

    if (m_bshowhomo)
    {
        if (m_beta)
            _listOrbitals2->clearSelection();

        disconnect(_listOrbitals, SIGNAL(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)),this,SLOT(OnOrbitalChange(QTreeWidgetItem*)));
        _listOrbitals->clearSelection();
        _listOrbitals->setCurrentItem(_listOrbitals->topLevelItem(m_render.OrbitalsData().Homo()-1));
        connect(_listOrbitals, SIGNAL(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)),this,SLOT(OnOrbitalChange(QTreeWidgetItem*)));

        m_render.ShowHomo();
        m_bonshow = true;
    }
    else
    {
        m_bonshow = false;
    }

    emit drawDensity(m_bonshow);

    if ((m_bonshow)&&(!_contoursGroupBox->isHidden()))
        m_plotspectrogram->FillSpectrogram();

}

void QOrbitalWidget::OnShowLumo()
{
    m_bshowlumo = !m_bshowlumo;
    m_bshowtotaldensity = false;
    m_bshowspindensity = false;
    m_bshowhomo = false;
    m_orbital = 0;//m_render.OrbitalsData().Lumo();
    m_betaorbital = 0;

    emit offtransitions(true);

    if (m_bshowlumo)
    {
        if (m_beta)
        {
            disconnect(_listOrbitals2, SIGNAL(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)),this,SLOT(OnOrbitalChange(QTreeWidgetItem*)));
            _listOrbitals2->clearSelection();
            _listOrbitals2->setCurrentItem(_listOrbitals2->topLevelItem(m_render.OrbitalsData().Lumo()-1));
            connect(_listOrbitals2, SIGNAL(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)),this,SLOT(OnOrbitalChange(QTreeWidgetItem*)));

            _listOrbitals->clearSelection();
        }
        else
        {
            disconnect(_listOrbitals, SIGNAL(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)),this,SLOT(OnOrbitalChange(QTreeWidgetItem*)));
            _listOrbitals->clearSelection();
            _listOrbitals->setCurrentItem(_listOrbitals->topLevelItem(m_render.OrbitalsData().Lumo()-1));
            connect(_listOrbitals, SIGNAL(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)),this,SLOT(OnOrbitalChange(QTreeWidgetItem*)));
        }

        m_render.ShowLumo();
        m_bonshow = true;
    }
    else
    {
        m_bonshow = false;
    }

    emit drawDensity(m_bonshow);

    if ((m_bonshow)&&(!_contoursGroupBox->isHidden()))
        m_plotspectrogram->FillSpectrogram();
}

void QOrbitalWidget::OnShowTotalDensity()
{
    m_bshowtotaldensity = !m_bshowtotaldensity;
    m_bshowspindensity = false;
    m_bshowhomo = false;
    m_bshowlumo = false;
    m_orbital = 0;
    m_betaorbital = 0;

    emit offtransitions(true);

    if (m_bshowtotaldensity)
    {
        if (m_beta)
            _listOrbitals2->clearSelection();
        _listOrbitals->clearSelection();

        m_render.ShowTotalDensity();
        m_bonshow = true;
    }
    else
        m_bonshow = false;

    emit drawDensity(m_bonshow);

    if ((m_bonshow)&&(!_contoursGroupBox->isHidden()))
        m_plotspectrogram->FillSpectrogram();
}

void QOrbitalWidget::OnShowSpinDensity()
{
    m_bshowspindensity = !m_bshowspindensity;
    m_bshowhomo = false;
    m_bshowlumo = false;
    m_bshowtotaldensity = false;
    m_orbital = 0;
    m_betaorbital = 0;

    emit offtransitions(true);

    if (m_bshowspindensity)
    {
        _listOrbitals2->clearSelection();
        _listOrbitals->clearSelection();

        m_render.ShowSpinDensity();
        m_bonshow = true;
    }
    else
        m_bonshow = false;

    emit drawDensity(m_bonshow);

    if ((m_bonshow)&&(!_contoursGroupBox->isHidden()))
        m_plotspectrogram->FillSpectrogram();
}

void QOrbitalWidget::OnOrbitalChange(QTreeWidgetItem* item)
{
    QTime timer;
    timer.start();
    m_bshowtotaldensity = false;
    m_bshowspindensity = false;
    m_bshowhomo = false;
    m_bshowlumo = false;
    m_orbital = item->text(0).toInt();
    m_betaorbital = 0;

    emit offtransitions(true);

    if (m_orbital>0)
    {
        qDebug() << "orbital=" << m_orbital;
        qDebug()  << "beta=" << m_beta;
        if (m_beta)
            _listOrbitals2->clearSelection();

        disconnect(_listOrbitals, SIGNAL(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)),this,SLOT(OnOrbitalChange(QTreeWidgetItem*)));
        _listOrbitals->clearSelection();
        _listOrbitals->setCurrentItem(_listOrbitals->topLevelItem(m_orbital-1));
        connect(_listOrbitals, SIGNAL(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)),this,SLOT(OnOrbitalChange(QTreeWidgetItem*)));

        m_render.ShowMolecularOrbital(m_orbital-1);
        m_bonshow = true;
    }
    else
        m_bonshow = false;

    emit drawDensity(m_bonshow);

    std::cout << "time needed is" << timer.elapsed()  << std::endl;
    if ((m_bonshow)&&(!_contoursGroupBox->isHidden()))
        m_plotspectrogram->FillSpectrogram();
}


void QOrbitalWidget::OnBetaOrbitalChange(QTreeWidgetItem* item)
{
    m_bshowtotaldensity = false;
    m_bshowspindensity = false;
    m_bshowhomo = false;
    m_bshowlumo = false;
    m_orbital = 0;
    m_betaorbital = item->text(0).toInt();

    _listOrbitals->clearSelection();

    emit offtransitions(true);

    qDebug() << "orbital beta=" << m_betaorbital;
    if (m_betaorbital > 0)
    {
        disconnect(_listOrbitals2, SIGNAL(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)),this,SLOT(OnOrbitalChange(QTreeWidgetItem*)));
        _listOrbitals2->clearSelection();
        _listOrbitals2->setCurrentItem(_listOrbitals2->topLevelItem(m_betaorbital-1));
        connect(_listOrbitals2, SIGNAL(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)),this,SLOT(OnOrbitalChange(QTreeWidgetItem*)));

        m_render.ShowBetaMolecularOrbital(m_betaorbital-1);
        m_bonshow = true;
    }
    else
        m_bonshow = false;

    emit drawDensity(m_bonshow);

    if ((m_bonshow)&&(!_contoursGroupBox->isHidden()))
        m_plotspectrogram->FillSpectrogram();
}

void QOrbitalWidget::OnShowTransitionChange(int item)
{
    m_bshowtotaldensity = false;
    m_bshowspindensity = false;
    m_bshowhomo = false;
    m_bshowlumo = false;
    m_orbital = 0;
    m_betaorbital = 0;

    _listOrbitals->clearSelection();
    if (m_beta)
        _listOrbitals2->clearSelection();

    if (item>=0)
    {
        m_render.ShowTransitionChange(item);
        m_bonshow = true;
    }
    else
        m_bonshow = false;

    emit drawDensity(m_bonshow);

    if ((m_bonshow)&&(!_contoursGroupBox->isHidden()))
        m_plotspectrogram->FillSpectrogram();

}

void QOrbitalWidget::OnShowDensityChange(int item)
{
    m_bshowtotaldensity = false;
    m_bshowspindensity = false;
    m_bshowhomo = false;
    m_bshowlumo = false;
    m_orbital = 0;
    m_betaorbital = 0;

    _listOrbitals->clearSelection();
    if (m_beta)
        _listOrbitals2->clearSelection();

    if (item>=0)
    {
        m_render.ShowDensityChange(item);
        m_bonshow = true;
    }
    else
        m_bonshow = false;

    emit drawDensity(m_bonshow);

    if ((m_bonshow)&&(!_contoursGroupBox->isHidden()))
        m_plotspectrogram->FillSpectrogram();
}


void QOrbitalWidget::OnShowContoursFrame()
{
    qDebug() << m_bshowcontours << endl;
    m_bshowcontours = !m_bshowcontours;
    qDebug() << m_bshowcontours << endl;

    if (m_bshowcontours)
    {

        if (!m_render.DensityData().DensityMatrix().Empty())
        {
            m_plotspectrogram->FillSpectrogram();
            _contoursGroupBox->show();
        }
    }
    else
        _contoursGroupBox->hide();
}


void QOrbitalWidget::OnChangeAxisPlane(bool b)
{
    if (_checkBoxXY->isChecked())
    {
        _spinBoxXY->setEnabled(true);
        _spinBoxXZ->setEnabled(false);
        _spinBoxYZ->setEnabled(false);
        m_plotspectrogram->SetAxisPlane(QPlotSpectrogram::XY);
        m_plotspectrogram->SetAxisValue(_spinBoxXY->value());
    }
    if (_checkBoxXZ->isChecked())
    {
        _spinBoxXZ->setEnabled(true);
        _spinBoxXY->setEnabled(false);
        _spinBoxYZ->setEnabled(false);
        m_plotspectrogram->SetAxisPlane(QPlotSpectrogram::XZ);
        m_plotspectrogram->SetAxisValue(_spinBoxXZ->value());
    }
    if (_checkBoxYZ->isChecked())
    {
        _spinBoxYZ->setEnabled(true);
        _spinBoxXZ->setEnabled(false);
        _spinBoxXY->setEnabled(false);
        m_plotspectrogram->SetAxisPlane(QPlotSpectrogram::YZ);
        m_plotspectrogram->SetAxisValue(_spinBoxYZ->value());
    }
    m_plotspectrogram->FillSpectrogram();

}


void QOrbitalWidget::OnChangeAxisValue(double v)
{
    m_plotspectrogram->SetAxisValue(v);
    m_plotspectrogram->FillSpectrogram();

}

void QOrbitalWidget::OffButtons(bool b)
{
    if (b)
    {
        m_bonshow = false;
        m_bshowtotaldensity = false;
        m_bshowspindensity = false;
        m_bshowhomo = false;
        m_bshowlumo = false;
        m_orbital = 0;
        m_betaorbital = 0;
        _listOrbitals->clearSelection();
        if (m_beta)
            _listOrbitals2->clearSelection();

        emit drawDensity(m_bonshow);
    }
}

void QOrbitalWidget::HideAllButtons()
{
    _labelIsovalue->hide();
    _sliderIsovalue->hide();
    _labelTransparence->hide();
    _sliderTransparence->hide();
    _labelGridResolution->hide();
    _sliderGridResolution->hide();
    _labelThreshold->hide();
    _thresholdSpinBox->hide();
    _labelOrbitals->hide();
    _labelOrbitals2->hide();
    _listOrbitals->hide();
    _listOrbitals2->hide();
    _buttonHomo->hide();
    _buttonLumo->hide();
    _buttonTotalDensity->hide();
    _buttonSpinDensity->hide();
    _buttonShowContours->hide();
    _contoursGroupBox->hide();

    m_bonshow = false;
    m_bshowtotaldensity = false;
    m_bshowspindensity = false;
    m_bshowhomo = false;
    m_bshowlumo = false;
    m_orbital = 0;
    m_betaorbital = 0;

    _listOrbitals->clearSelection();
    if (m_beta)
        _listOrbitals2->clearSelection();
}

void QOrbitalWidget::OnSetFrame(size_t frame)
{
    Q_ASSERT(m_world);
    if ( frame < 0 || frame >= m_world->Molecules().back().Frames().size() ) return;
    /*QMessageBox::information(nullptr,"title","frame: "+QString::number(frame+1)
                             +"\n"+"orbital: "+QString::number(m_orbital));*/
    kryomol::Frame& fr=m_world->Molecules().back().Frames()[frame];
    if ( fr.HasOrbitals() )
    {
        this->ListOrbitals();
        m_render = kryomol::RenderOrbitals(m_world->Molecules().back().Frames()[frame]);


        if (m_orbital > 0)
        {
            if (m_beta)
                _listOrbitals2->clearSelection();

            disconnect(_listOrbitals, SIGNAL(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)),this,SLOT(OnOrbitalChange(QTreeWidgetItem*)));
            _listOrbitals->clearSelection();
            _listOrbitals->setCurrentItem(_listOrbitals->topLevelItem(m_orbital-1));
            connect(_listOrbitals, SIGNAL(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)),this,SLOT(OnOrbitalChange(QTreeWidgetItem*)));

            m_render.ShowMolecularOrbital(m_orbital-1);
            m_bonshow = true;
        } else m_bonshow=false;
        emit drawDensity(m_bonshow);
    }
}

void QOrbitalWidget::OnDrawDensity(bool b)
{
    if (b)
    {
        m_world->CurrentMolecule()->CurrentFrame().SetPositiveDensity(this->GetDensity().PositiveRenderDensityVector());
        m_world->CurrentMolecule()->CurrentFrame().SetNegativeDensity(this->GetDensity().NegativeRenderDensityVector());

    }
    m_world->OnShowDensity(b);
}
