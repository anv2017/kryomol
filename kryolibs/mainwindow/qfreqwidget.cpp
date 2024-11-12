/*****************************************************************************************
                            qfreqwidget.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <iostream>
#include <sstream>

#include <QTableView>
#include <QHeaderView>
#include <qpushbutton.h>
#include <qfiledialog.h>
#include <qlineedit.h>
#include <qmessagebox.h>
#include <qvalidator.h>
#include <qslider.h>
#include <qspinbox.h>
#include <qcombobox.h>
#include <QFrame>

#include "qfreqwidget.h"
#include "irspectrum.h"
#include "qlogtablewidget.h"
#include "qplotspectrum.h"

#include "qryomolapp.h"

#include "world.h"
#include "molecule.h"
#include "frame.h"

using namespace kryomol;

QFreqWidget::QFreqWidget(kryomol::World* world,const QString& file,QWidget* parent)
    : QWidget(parent), m_world(world)
{
    setupUi(this);
    if( file != QString::null)
        m_file=(std::string) file.toLatin1();
    else
        m_file="Unnamed";

    m_bshowspectrum=true;

    _lineWidthSpinBox->setValue(LineWidth());
    _shift->setValue(Shift());
    _numpoints->setValue(NPoints());

    QDoubleValidator* val= new QDoubleValidator(this);
    val->setBottom(0);
    val->setTop(1);
    _scaleLineEdit->setValidator(val);
    _scaleLineEdit->setText("1.00");

    m_freqtable= new QLogTableWidget(this);
    m_freqtable->setSortingEnabled(true);
    m_freqtable->setCopyHorizontalHeader(true);
    m_freqtable->verticalHeader()->hide();
    QVBoxLayout* tablelay=new QVBoxLayout(_freqFrame);
    tablelay->addWidget(m_freqtable);

    //QHeaderView* vheader= m_freqtable->verticalHeader();
    //connect(vheader,SIGNAL(sectionClicked(int)),this,SLOT(OnTableSelection(int)));
    //connect(m_freqtable->horizontalHeader(),SIGNAL(sectionClicked(int)),m_freqtable,SLOT( sortColumn(int)));
    connect(m_freqtable,SIGNAL(cellClicked(int,int)),this,SLOT(OnTableSelection(int)));
    connect(showSpecButton,SIGNAL(clicked()),this,SLOT(OnShowSpectrum()));
    connect(_exportIRButton,SIGNAL(clicked()),this,SLOT(OnWriteJCampDX()));
    connect(_copyButton,SIGNAL(clicked()),this,SLOT(OnCopy()));
    connect(_minLimitLineEdit,SIGNAL(returnPressed()),this,SLOT(OnChangeLimits()));
    connect(_maxLimitLineEdit,SIGNAL(returnPressed()),this,SLOT(OnChangeLimits()));
    connect(_resetLimitsButton,SIGNAL(clicked()),this,SLOT(OnResetLimits()));
    connect(_distortionSlider,SIGNAL(valueChanged(int)),this,SLOT(OnDistort(int)));
    connect(_resetDistortionButton,SIGNAL(clicked()),this,SLOT(OnResetDistortions()));
    connect(_lineWidthSpinBox,SIGNAL(valueChanged(double)),this,SLOT(OnSetLineWidth(double)));
    connect(_shift,SIGNAL(valueChanged(int)),this,SLOT(OnSetShift(int)));
    connect(_numpoints,SIGNAL(valueChanged(int)),this,SLOT(OnSetNPoints(int)));

    connect(_scaleLineEdit,SIGNAL(returnPressed()),this,SLOT(OnScaling()));
    connect(_spectrumTypeComboBox,SIGNAL(activated(int)),this,SLOT(OnSpectrumTypeChanged(int)));

    m_activemode=0;
    m_npoints = 8192;
}


QFreqWidget::~QFreqWidget()
{}

//Init table for conformer with index fidx
void QFreqWidget::InitTable(size_t fidx)
{

    if ( m_frequencysets.empty() )
    {
        std::cout << "Empty list of frequencies. Table with not be initialized";
        return;
    }

    QString sheader;
    if (GetType() == QPlotSpectrum::IR)
        sheader="IR Intensity (kM/mol)";
    if (GetType() == QPlotSpectrum::VCD)
        sheader="Rotatory Strength";
    if (GetType() == QPlotSpectrum::RAMAN)
        sheader="Activity";

    QStringList hheader;
    hheader << "Mode";
    hheader << "Freq(cm-1)";
    hheader << sheader;

    const std::vector<Frequency>& frequencies=m_frequencysets[fidx];

    std::vector<Frequency>::iterator it;

    m_freqtable->setRowCount(frequencies.size());
    m_freqtable->setColumnCount(3);
    m_distortframes.resize(frequencies.size());
    for(size_t i=0;i<frequencies.size();i++)
    {
        m_distortframes[i]=0;
    }

    //QStringList vheader;

    int i;
    for(size_t i=0;i<frequencies.size();++i)
    {
        auto& f=frequencies[i];
        //QString sf, si;
        //sf.setNum(it->x,'f',1);
        float yvalue;
        int precision;
        switch(GetType())
        {
        case QPlotSpectrum::VCD:
            //si.setNum(it->z,'f',3);
            yvalue=f.z;
            precision=3;
            break;
        case QPlotSpectrum::IR:
            //si.setNum(it->y,'f',1);
            yvalue=f.y;
            precision=1;
            break;
        case QPlotSpectrum::RAMAN:
            yvalue=f.w;
            precision=1;
            break;
        default:
            throw QString("Not handled QPlotSpectrum type");
            break;
        }

        // item=new QTableWidgetItem ( sf);
        // item->setFlags ( Qt::ItemIsEnabled | Qt::ItemIsSelectable );
        m_freqtable->setNum ( i,1,f.x,'f',1 );

        //item=new QTableWidgetItem ( si);
        //item->setFlags ( Qt::ItemIsEnabled | Qt::ItemIsSelectable );
        m_freqtable->setNum ( i,2,yvalue,'f',precision);

        //    QString str;
        //    str.sprintf("%d",i+1);
        //    vheader << str;

        m_freqtable->setNum(i,0,static_cast<double>(i+1));
        m_freqtable->item(i,0)->setBackground(QApplication::palette().button());

    }

    m_freqtable->setHorizontalHeaderLabels ( hheader );
    //m_freqtable->setVerticalHeaderLabels ( vheader );

    m_freqtable->horizontalHeader()->setSortIndicator(0,Qt::AscendingOrder);
    m_freqtable->horizontalHeader()->setSortIndicatorShown(true);

    m_freqtable->resizeColumnsToContents();

    //initialize the irspectrum object
    SetAutoLimits();

    QString max,min;
    max.setNum(Max());
    min.setNum(Min());
    _minLimitLineEdit->setText(min);
    _maxLimitLineEdit->setText(max);

    CalculateSpectrum();
    //emit Type(GetType()); //to stablish the baseline
    _numpoints->setValue(m_npoints);
}


void QFreqWidget::OnWriteJCampDX()
{
    //propose a name for the file
    QString sfile(m_file.c_str());
    QString sfilename=sfile.section('.',0,-2)+".jdx"; //change out extension by jdx
    QString file= QFileDialog::getSaveFileName(this,"Write JCAMP-DX spectrum",sfilename,"JCAMP-DX (*.jdx)");//(sfilename,"JCAMP-DX (*.jdx)",this,"KryoMoll","Write JCAMP-DX spectrum" );
    if( file != QString::null)
    {
        m_file=(std::string) file.toLatin1();
        WriteJCampDX();
    }
}


void QFreqWidget::OnCopy()
{
    CopyData();
}


void QFreqWidget::OnChangeLimits()
{
    float min =_minLimitLineEdit->text().toFloat();
    float max =_maxLimitLineEdit->text().toFloat();
    if( max < min)
    {
        QMessageBox mb;
        mb.setText("Sorry, limits are inverted");
        mb.exec();
        return;
    }
    SetLimits(max,min);
    CalculateSpectrum();
    emit data(GetData(),Max(),Min(),Shift());
}


void QFreqWidget::OnResetLimits()
{
    SetAutoLimits();
    //write the limits
    QString max,min;
    max.setNum(Max());
    min.setNum(Min());
    _minLimitLineEdit->setText(min);
    _maxLimitLineEdit->setText(max);
    CalculateSpectrum();
    emit data(GetData(),Max(),Min(),Shift());
}


void QFreqWidget::OnSetLimits(float max,float min)
{
    QString smax,smin;
    smax.setNum(max);
    smin.setNum(min);
    _minLimitLineEdit->setText(smin);
    _maxLimitLineEdit->setText(smax);
    SetLimits(max,min);
    CalculateSpectrum();
    emit data(GetData(),Max(),Min(),Shift());
}


void QFreqWidget::OnSetLineWidth(double lw)
{
    SetLineWidth(lw);
    CalculateSpectrum();
    emit data(GetData(),Max(),Min(),Shift());
}


void QFreqWidget::OnSetShift(int lw)
{
    SetShift(lw);
    CalculateSpectrum();
    emit data(GetData(),Max(),Min(),Shift());
}

void QFreqWidget::OnSetNPoints(int lw)
{
    if (lw>m_npoints)
    {
        m_npoints = m_npoints*2;
        _numpoints->setValue(m_npoints);
    }
    else
    {
        if (lw<m_npoints)
        {
            m_npoints = m_npoints/2;
            _numpoints->setValue(m_npoints);
        }
    }
    SetNPoints(m_npoints);
    CalculateSpectrum();
    emit data(GetData(),Max(),Min(),Shift());
}

void QFreqWidget::OnShowSpectrum()
{
    emit showspectrum(m_bshowspectrum);
    emit Type(GetType()); //to stablish the baseline
    emit data(GetData(),Max(),Min(),Shift());

    m_bshowspectrum=!m_bshowspectrum;

    QString buttontext;
    if(!m_bshowspectrum)
        buttontext=tr("Hide Spectrum");
    else
        buttontext=tr("Show Spectrum");

    showSpecButton->setText(buttontext);
}


void QFreqWidget::OnScaling()
{
    float scale=_scaleLineEdit->text().toFloat();

    std::vector<Frequency>::iterator it;
    IRSpectrum::SetFrequencies(m_frequencysets,scale);

    /*int i;
  for(it=m_frequencies.begin(),i=0;it!=m_frequencies.end();it++,i++)
  {
   // QString sf, si;
   // sf.setNum(it->x,'f',1);
   // si.setNum(it->y,'f',1);

   // QTableWidgetItem* item=new QTableWidgetItem ( sf);
    //item->setFlags ( Qt::ItemIsEnabled | Qt::ItemIsSelectable );
    m_freqtable->setNum ( i,0,it->x,'f',1 );

    //item=new QTableWidgetItem ( si);
    //item->setFlags ( Qt::ItemIsEnabled | Qt::ItemIsSelectable );
    m_freqtable->setNum ( i,1,it->y,'f',1 );
  }*/
    InitTable(m_world->CurrentMolecule()->CurrentFrameIndex());

    //initialize the irspectrum object
    SetAutoLimits();

    QString max,min;
    max.setNum(Max());
    min.setNum(Min());
    _minLimitLineEdit->setText(min);
    _maxLimitLineEdit->setText(max);

    CalculateSpectrum();
    emit data(GetData(),Max(),Min(),Shift());
}


void QFreqWidget::InitFrequencies()
{

    for(const auto& f : m_world->CurrentMolecule()->Frames() )
    {
        m_frequencysets.push_back(f.GetFrequencies());
    }

    //the above set will not be modified, scaled etc but the irspectrum set can be scaled, shifted, etc.
    IRSpectrum::SetFrequencies(m_frequencysets);
}


void QFreqWidget::OnTableSelection(int s)
{
    if ( s < 0 || s > m_freqtable->rowCount()-1 ) return;
    int signal = m_freqtable->item(s,0)->text().toInt();
    OnAnimate(signal-1);
}


void QFreqWidget::OnAnimate(int m)
{
    m_activemode=m;
    _distortionSlider->blockSignals(true);
    _distortionSlider->setValue(m_distortframes[m_activemode]);
    _distortionSlider->blockSignals(false);
    emit mode(m,m_distortframes[m]);
}


void QFreqWidget::OnDisableControls(bool b)
{
    _distortionSlider->setEnabled(!b);
    _resetDistortionButton->setEnabled(!b);
}


void QFreqWidget::OnDistort(int d)
{
    m_distortframes[m_activemode]=_distortionSlider->value();
    emit distort(d);
}

void QFreqWidget::OnResetDistortions()
{
    emit reset();
    for(size_t i=0;i<m_distortframes.size();i++)
    {
        m_distortframes[i]=0;
    }
    _distortionSlider->blockSignals(true);
    _distortionSlider->setValue(0);
    _distortionSlider->blockSignals(false);
}


void QFreqWidget::OnSpectrumTypeChanged(int t)
{
    switch(t)
    {
    case 0:
       SetType(QPlotSpectrum::IR);
        break;
    case 1:
       SetType(QPlotSpectrum::VCD);
        break;
    case 2:
        SetType(QPlotSpectrum::RAMAN);
        break;
    default:
        throw QString("Invalid type");
    }

    InitTable(m_world->CurrentMolecule()->CurrentFrameIndex());
    emit data(GetData(),Max(),Min(),Shift());
    emit Type(GetType()); //to stablish the baseline
}


