/*****************************************************************************************
                            quvwidget.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <iostream>

#include <qtablewidget.h>
#include <qlineedit.h>
#include <qspinbox.h>
#include <qpushbutton.h>
#include <qcombobox.h>
#include <qmessagebox.h>
#include <qfiledialog.h>
#include <QHeaderView>

#include "quvwidget.h"
#include "qlogtable.h"
#include "qlogtablewidget.h"

#include "frequency.h"
#include "world.h"
#include "molecule.h"
#include "frame.h"

QUVWidget::QUVWidget(const kryomol::World* world, const QString& file, QWidget* parent)
    : QWidget(parent) , m_world(world)
{
    setupUi(this);
    if( file != QString::null)
        m_file=(std::string) file.toLatin1();
    else
        m_file="Unnamed";

    m_uvTable= new QLogTableWidget(this);
    m_uvTable->setSortingEnabled(true);
    m_uvTable->setCopyHorizontalHeader(true);
    m_uvTable->verticalHeader()->hide();
    m_tablelay=new QVBoxLayout(_uvFrame);
    //tablelay->addWidget(m_uvTable);



    m_bshowspectrum=true;
    connect(showSpecButton,SIGNAL(clicked()),this,SLOT(OnShowSpectrum()));
    connect(_spectrumTypeComboBox,SIGNAL(currentIndexChanged(int)),this,SLOT(OnSpectrumTypeChanged(int)));

    _lineWidthSpinBox->setValue(LineWidth());
    _shift->setValue(Shift());
    _numpoints->setValue(NPoints());
    connect(_lineWidthSpinBox,SIGNAL(valueChanged(double)),this,SLOT(OnSetLineWidth(double)));
    connect(_resetLimitsButton,SIGNAL(clicked()),this,SLOT(OnResetLimits()));
    connect(_exportIRButton,SIGNAL(clicked()),this,SLOT(OnWriteJCampDX()));
    connect(_copyButton,SIGNAL(clicked()),this,SLOT(OnCopy()));
    connect(_shift,SIGNAL(valueChanged(int)),this,SLOT(OnSetShift(int)));
    connect(_numpoints,SIGNAL(valueChanged(int)),this,SLOT(OnSetNPoints(int)));

    m_bshowcoefficients = false;
    m_bshowtransitions = false;
    m_bshowdensities = false;
    m_bshowelectricdipole = false;
    m_bshowmagneticdipole = false;
    m_bshowvelocitydipole = false;
    connect(_leftLimitLineEdit,SIGNAL(returnPressed()),this,SLOT(OnChangeLimits()));
    connect(_rightLimitLineEdit,SIGNAL(returnPressed()),this,SLOT(OnChangeLimits()));
    connect(_electricDipoleCheckBox,SIGNAL(toggled(bool)),this,SLOT(OnShowElectricDipole(bool)));
    connect(_magneticDipoleCheckBox,SIGNAL(toggled(bool)),this,SLOT(OnShowMagneticDipole(bool)));
    connect(_velocityDipoleCheckBox,SIGNAL(toggled(bool)),this,SLOT(OnShowVelocityDipole(bool)));
    connect(_densityChangesCheckBox,SIGNAL(toggled(bool)),this,SLOT(OnShowDensityChanges(bool)));
    connect(_transitionChangesCheckBox,SIGNAL(toggled(bool)),this,SLOT(OnShowTransitionChanges(bool)));
    connect(_transitionCoefficientsCheckBox,SIGNAL(toggled(bool)),this,SLOT(OnShowTransitionCoefficients(bool)));

    m_bsolventshift=false;
    connect(_solventShiftCheckBox,SIGNAL(toggled(bool)),this,SLOT(OnSolventShift(bool )));

    _densityChangesCheckBox->setDisabled(true);
    _transitionChangesCheckBox->setDisabled(true);
    _transitionCoefficientsCheckBox->setDisabled(true);

    _transitionTable->hide();

    _enantiomerButton->hide();
    _enantiomerButton->setChecked(Enantiomer());
    _ECDFormalismComboBox->addItem("Length");
    _ECDFormalismComboBox->addItem("Velocity");
    connect(_ECDFormalismComboBox,SIGNAL(currentIndexChanged(int)),this,SLOT(OnFormalismChanged(int)));
    _ECDFormalismComboBox->hide();
    connect(_enantiomerButton,SIGNAL(toggled(bool)),this,SLOT(OnCalculateEnantiomer(bool )));

}




QUVWidget::~QUVWidget()
{}

void QUVWidget::OnShowSpectrum()
{
    emit showspectrum(m_bshowspectrum);
    emit type(GetType()); //to stablish the baseline
    emit data(GetData(),GetTotalData(),Max(),Min(),Shift());

    m_bshowspectrum=!m_bshowspectrum;

    QString buttontext;
    if(!m_bshowspectrum)
    {
        buttontext=tr("Hide Spectrum");
    }
    else buttontext=tr("Show Spectrum");

    showSpecButton->setText(buttontext);
}


void QUVWidget::InitTable(size_t fidx)
{
    const std::vector<Spectralline>& lines=m_linesets.at(fidx);
    delete m_uvTable;
    m_uvTable  = new QLogTableWidget(this);
    m_uvTable->setSortingEnabled(true);
    m_uvTable->setCopyHorizontalHeader(true);
    m_uvTable->verticalHeader()->hide();
    m_tablelay->addWidget(m_uvTable);
    m_uvTable->setRowCount(lines.size());
    connect(m_uvTable,SIGNAL(cellClicked(int,int)),this,SLOT(OnTransitionSelected(int )));
    QStringList hheader;
    if (GetType() == QPlotSpectrum::UV)
    {
        hheader << "Transition" << "Lambda(nm)" << "Oscillator Strenght";
        m_uvTable->setColumnCount(3);
    }
    if (GetType() == QPlotSpectrum::ECD)
    {
        hheader << "Transition" << "Lambda(nm)" << "Rotatory strength (length )" << "Rotatory strenght (velocity) ";
        m_uvTable->setColumnCount(4);
    }

    for(int i=0;i<(int)lines.size();++i)
    {
        const Spectralline& line= lines.at(i);
        float yvalue1, yvalue2, yvalue3;
        if ( GetType() == QPlotSpectrum::UV )
        {
            yvalue1=line.y0;
        }
        if (GetType() == QPlotSpectrum::ECD)
        {
            yvalue1=line.RotatoryStrengthLength();
            yvalue2=line.RotatoryStrengthVelocity();
        }

        m_uvTable->setNum(i,0,i+1);
        float x=line.x;
        if ( m_bsolventshift )
        {
            x-=line.SolventShift();
        }
        m_uvTable->setNum(i,1,x);
        m_uvTable->setNum ( i,2,yvalue1);

        if ( GetType() == QPlotSpectrum::ECD  )
            m_uvTable->setNum( i,3, yvalue2 );

        m_uvTable->item(i,0)->setBackground(QApplication::palette().button());

    }

    m_uvTable->setHorizontalHeaderLabels ( hheader );
    m_uvTable->horizontalHeader()->setSortIndicator(0,Qt::AscendingOrder);
    m_uvTable->horizontalHeader()->setSortIndicatorShown(true);
    m_uvTable->resizeColumnsToContents();


    if ((m_bshowcoefficients) && (m_uvTable->currentRow()>=0))
    {

        int transitionindex=m_uvTable->item(m_uvTable->currentRow(),0)->text().toInt()-1;
        InitTransitionTable(transitionindex);

    }

    //initialize the irspectrum object

    SetAutoLimits();
    //write the limits
    QString max,min;
    max.setNum(Max());
    min.setNum(Min());
    _leftLimitLineEdit->setText(max);
    _rightLimitLineEdit->setText(min);
    CalculateSpectrum();

}

void QUVWidget::InitTransitionTable(int transition)
{
    size_t fidx=this->m_world->CurrentMolecule()->CurrentFrameIndex();
    _transitionTable->clear();

    _transitionTable->setColumnCount(3);

    QStringList hheader;
    hheader << "Coeff" << "MO(i)" << "MO(j)";

    QStringList vheader;

    const std::vector<kryomol::TransitionChange>& transitions = m_transitions.at(fidx).at(transition);
    _transitionTable->setRowCount(transitions.size());

    QTableWidgetItem* item;
    for(size_t i=0;i<transitions.size();++i)
    {
        QString coef, orb_i, orb_j;
        coef.setNum(transitions.at(i).Coefficient(),'f',3);
        item = new QTableWidgetItem (coef);
        item->setFlags ( Qt::ItemIsEnabled | Qt::ItemIsSelectable );
        _transitionTable->setItem ( i,0,item );

        if (m_beta)
        {
            orb_i = QString(transitions.at(i).OrbitalSI().c_str());
            orb_j = QString(transitions.at(i).OrbitalSJ().c_str());
        }
        else
        {
            orb_i.setNum(transitions.at(i).OrbitalI());
            orb_j.setNum(transitions.at(i).OrbitalJ());
        }
        item = new QTableWidgetItem (orb_i);
        item->setFlags ( Qt::ItemIsEnabled | Qt::ItemIsSelectable );
        _transitionTable->setItem ( i,1,item );

        item = new QTableWidgetItem (orb_j);
        item->setFlags ( Qt::ItemIsEnabled | Qt::ItemIsSelectable );
        _transitionTable->setItem ( i,2,item );

        QString str;
        str.sprintf("%d",i+1);
        vheader << str;
    }


    _transitionTable->setHorizontalHeaderLabels ( hheader );
    _transitionTable->setVerticalHeaderLabels ( vheader );
    _transitionTable->horizontalHeader()->setSortIndicator(0,Qt::AscendingOrder);
    _transitionTable->horizontalHeader()->setSortIndicatorShown(true);
    _transitionTable->resizeColumnsToContents();


}

void QUVWidget::OnTransitionSelected(int row)
{
    int tindex=m_uvTable->item(row,0)->text().toInt() - 1 ;
    if (m_bshowcoefficients)
        InitTransitionTable(tindex);

    emit transitionselected( tindex );

    emit setactivetransition(tindex);
    if (m_bshowelectricdipole)
        emit(showelectricdipole(true));
    if (m_bshowmagneticdipole)
        emit(showmagneticdipole(true));
    if (m_bshowvelocitydipole)
        emit(showvelocitydipole(true));

    if (m_bshowtransitions)
        emit showtransition( tindex );

    if (m_bshowdensities)
        emit showdensities( tindex );
}


void QUVWidget::OnSpectrumTypeChanged(int t)
{

    if (t==0)
    {
        SetType(QPlotSpectrum::UV);
        _enantiomerButton->hide();
        _ECDFormalismComboBox->hide();
    }
    else
    {
        SetType(QPlotSpectrum::ECD);
        _enantiomerButton->show();
        _ECDFormalismComboBox->show();
    }
    InitTable(this->m_world->CurrentMolecule()->CurrentFrameIndex());
    emit data(GetData(),GetTotalData(),Max(),Min(),Shift());
    emit type(GetType());
}

void QUVWidget::OnSetLineWidth(double lw)
{
    SetLineWidth(lw);
    CalculateSpectrum();
    emit data(GetData(),GetTotalData(),Max(),Min(),Shift());
}

void QUVWidget::OnSetShift(int lw)
{
    SetShift(lw);
    CalculateSpectrum();
    emit data(GetData(),GetTotalData(),Max(),Min(),Shift());
}

void QUVWidget::OnSetNPoints(int lw)
{
    SetNPoints(lw);
    CalculateSpectrum();
    emit data(GetData(),GetTotalData(),Max(),Min(),Shift());
}


void QUVWidget::OnResetLimits()
{
    SetAutoLimits();
    //write the limits
    QString max,min;
    max.setNum(Max());
    min.setNum(Min());
    _leftLimitLineEdit->setText(max);
    _rightLimitLineEdit->setText(min);
    CalculateSpectrum();
    this->UpdatePlot();

}

void QUVWidget::UpdatePlot()
{
    emit data(GetData(),GetTotalData(),Max(),Min(),Shift());
}

void QUVWidget::OnChangeLimits()
{
    float min =_rightLimitLineEdit->text().toFloat();
    float max =_leftLimitLineEdit->text().toFloat();
    if( max < min)
    {
        QMessageBox mb;
        mb.setText("Sorry, limits are inverted");
        mb.exec();
        return;
    }
    SetLimits(max,min);
    CalculateSpectrum();
    UpdatePlot();
}

void QUVWidget::OnSetLimits(float max,float min)
{
    QString smax,smin;
    smax.setNum(max);
    smin.setNum(min);
    _leftLimitLineEdit->setText(smax);
    _rightLimitLineEdit->setText(smin);
    SetLimits(max,min);
    CalculateSpectrum();
    UpdatePlot();
}

void QUVWidget::OnWriteJCampDX()
{

#ifdef USE_KDE
#else
    //propose a name for the file
    QString sfile(m_file.c_str());
    QString sfilename=sfile.section('.',0,-2)+".jdx"; //change out extension by jdx
    QString file= QFileDialog::getSaveFileName(this,"Write JCAMP-DX spectrum",sfilename,"JCAMP-DX (*.jdx)");//this,"QryoMol", "Write JCAMP-DX spectrum" );
    if( file != QString::null)
    {
        m_file=(std::string) file.toLatin1();
        WriteJCampDX();
    }
#endif

}

void QUVWidget::OnCopy()
{
    CopyData();
}

void QUVWidget::OnShowDensityChanges(bool b)
{
    m_bshowdensities = !m_bshowdensities;
    if  (!m_bshowdensities)
        emit offshowtransitions(true);
    else
    {
        m_bshowtransitions = false;
        disconnect(_transitionChangesCheckBox,SIGNAL(toggled(bool)),this,SLOT(OnShowTransitionChanges(bool)));
        _transitionChangesCheckBox->setChecked(false);
        connect(_transitionChangesCheckBox,SIGNAL(toggled(bool)),this,SLOT(OnShowTransitionChanges(bool)));

        if (!m_uvTable->selectedItems().empty())
        {
            int transitionindex=m_uvTable->item(m_uvTable->currentRow(),0)->text().toInt()-1;
            emit showdensities(transitionindex);
        }
        else
            emit showdensities(-1);
    }
}

void QUVWidget::OnShowTransitionChanges(bool b)
{
    m_bshowtransitions = !m_bshowtransitions;
    if  (!m_bshowtransitions)
        emit offshowtransitions(true);
    else
    {
        m_bshowdensities = false;
        disconnect(_densityChangesCheckBox,SIGNAL(toggled(bool)),this,SLOT(OnShowDensityChanges(bool)));
        _densityChangesCheckBox->setChecked(false);
        connect(_densityChangesCheckBox,SIGNAL(toggled(bool)),this,SLOT(OnShowDensityChanges(bool)));

        if (!m_uvTable->selectedItems().empty())
        {
            int transitionindex=m_uvTable->item(m_uvTable->currentRow(),0)->text().toInt()-1;
            emit showtransition(transitionindex);
        }
        else
            emit showtransition(-1);
    }
}


void QUVWidget::OffShowTransitionChanges(bool b)
{
    m_bshowdensities = false;
    m_bshowtransitions = false;
    m_bshowcoefficients = false;

    _transitionTable->hide();

    emit showelectricdipole(false);
    emit showmagneticdipole(false);
    emit showvelocitydipole(false);

    disconnect(_electricDipoleCheckBox,SIGNAL(toggled(bool)),this,SIGNAL(showelectricdipole(bool)));
    disconnect(_magneticDipoleCheckBox,SIGNAL(toggled(bool)),this,SIGNAL(showmagneticdipole(bool)));
    disconnect(_velocityDipoleCheckBox,SIGNAL(toggled(bool)),this,SIGNAL(showvelocitydipole(bool)));
    disconnect(_densityChangesCheckBox,SIGNAL(toggled(bool)),this,SLOT(OnShowDensityChanges(bool)));
    disconnect(_transitionChangesCheckBox,SIGNAL(toggled(bool)),this,SLOT(OnShowTransitionChanges(bool)));
    disconnect(_transitionCoefficientsCheckBox,SIGNAL(toggled(bool)),this,SLOT(OnShowTransitionCoefficients(bool)));
    disconnect(m_uvTable,SIGNAL(cellClicked(int,int)),this,SLOT(OnTransitionSelected(int )));
    _densityChangesCheckBox->setChecked(false);
    _transitionChangesCheckBox->setChecked(false);
    _transitionCoefficientsCheckBox->setChecked(false);
    _electricDipoleCheckBox->setChecked(false);
    _magneticDipoleCheckBox->setChecked(false);
    _velocityDipoleCheckBox->setChecked(false);
    m_uvTable->clearSelection();
    connect(_electricDipoleCheckBox,SIGNAL(toggled(bool)),this,SIGNAL(showelectricdipole(bool)));
    connect(_magneticDipoleCheckBox,SIGNAL(toggled(bool)),this,SIGNAL(showmagneticdipole(bool)));
    connect(_velocityDipoleCheckBox,SIGNAL(toggled(bool)),this,SIGNAL(showvelocitydipole(bool)));
    connect(_densityChangesCheckBox,SIGNAL(toggled(bool)),this,SLOT(OnShowDensityChanges(bool)));
    connect(_transitionChangesCheckBox,SIGNAL(toggled(bool)),this,SLOT(OnShowTransitionChanges(bool)));
    connect(_transitionCoefficientsCheckBox,SIGNAL(toggled(bool)),this,SLOT(OnShowTransitionCoefficients(bool)));
    connect(m_uvTable,SIGNAL(cellClicked(int,int)),this,SLOT(OnTransitionSelected(int )));


}

void QUVWidget::OnSolventShift(bool b)
{
    m_bsolventshift=b;
    InitTable(this->m_world->CurrentMolecule()->CurrentFrameIndex());
    SubstractSolventShift(b);
    CalculateSpectrum();
    UpdatePlot();

}

void QUVWidget::SetCheckableTransitionChanges()
{
    _densityChangesCheckBox->setEnabled(true);
    _transitionChangesCheckBox->setEnabled(true);
    _transitionCoefficientsCheckBox->setEnabled(true);
    _densityChangesCheckBox->setCheckable(true);
    _transitionChangesCheckBox->setCheckable(true);
    _transitionCoefficientsCheckBox->setCheckable(true);
}

void QUVWidget::SetCheckableTransitionCoefficients()
{
    _transitionCoefficientsCheckBox->setEnabled(true);
    _transitionCoefficientsCheckBox->setCheckable(true);
}

void QUVWidget::OnShowTransitionCoefficients(bool b)
{
    m_bshowcoefficients = !m_bshowcoefficients;

    if (m_bshowcoefficients)
    {
        int transitionindex=m_uvTable->item(m_uvTable->currentRow(),0)->text().toInt()-1;
        InitTransitionTable(transitionindex);
        _transitionTable->show();
    }
    else
        _transitionTable->hide();
}

void QUVWidget::OnShowElectricDipole(bool b)
{
    m_bshowelectricdipole = !m_bshowelectricdipole;

    emit (showelectricdipole(m_bshowelectricdipole));
}

void QUVWidget::OnShowMagneticDipole(bool b)
{
    m_bshowmagneticdipole = !m_bshowmagneticdipole;

    emit (showmagneticdipole(m_bshowmagneticdipole));
}

void QUVWidget::OnShowVelocityDipole(bool b)
{
    m_bshowvelocitydipole = !m_bshowvelocitydipole;

    emit (showvelocitydipole(m_bshowvelocitydipole));
}

void QUVWidget::OnCalculateEnantiomer(bool b)
{
    SetEnantiomer(b);
    CalculateSpectrum();
    UpdatePlot();
}

void QUVWidget::OnFormalismChanged(int index)
{
    SetFormalism( static_cast<formalism>(index));
    CalculateSpectrum();
    UpdatePlot();
}

