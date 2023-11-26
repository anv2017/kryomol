/*****************************************************************************************
                            qmeasurewidget.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <iostream>
#include <sstream>

#include "qmeasurewidget.h"
#include "ui_qmeasurewidgetbase.h"

QMeasureWidget::QMeasureWidget(QWidget *parent) :
    QWidget(parent)
{
    setupUi ( this );

    Init();

}

QMeasureWidget::~QMeasureWidget()
{

}

void QMeasureWidget::Init()
{    
    connect(_clearAll, SIGNAL ( clicked() ), this, SLOT ( OnClearAll()) );
    connect(_listDistances, SIGNAL(itemActivated(QListWidgetItem*)),this,SLOT(OnDistanceChange(QListWidgetItem*)));
    connect(_listAngles, SIGNAL(itemActivated(QListWidgetItem*)),this,SLOT(OnAngleChange(QListWidgetItem*)));
    connect(_listDihedrals, SIGNAL(itemActivated(QListWidgetItem*)),this,SLOT(OnDihedralChange(QListWidgetItem*)));
    connect(_showDistances, SIGNAL (clicked()),this, SLOT (OnShowDistances()));

    m_bshowdistances=false;

    _clearAll->setWhatsThis(tr("Clear all the measurements calculated"));
    _showDistances->setWhatsThis(tr("Activate it for visualizing the distance vector in the visor"));
}

void QMeasureWidget::OnClearAll()
{
    _listDistances->clear();
    _listAngles->clear();
    _listDihedrals->clear();
    _showDistances->setCheckState(Qt::Unchecked);
    m_bshowdistances=false;

    emit clearAll();

}

void QMeasureWidget::OnWriteDistance(QString& str)
{
    new QListWidgetItem(str, _listDistances, 0);

}

void QMeasureWidget::OnWriteAngle(QString& str)
{
    new QListWidgetItem(str, _listAngles, 0);

}


void QMeasureWidget::OnWriteDihedral(QString& str)
{
    new QListWidgetItem(str, _listDihedrals, 0);

}

void QMeasureWidget::updateDistances(QStringList list)
{
    bool sel=false;
    int index;
    if (_listDistances->currentItem())
    {
        if (_listDistances->currentItem()->isSelected())
        {
            sel = true;
            index = _listDistances->currentIndex().row();
        }
    }
    _listDistances->clear();
    for (int i=0; i<list.count(); i++)
    {
        _listDistances->addItem(list.at(i));
    }
    if (sel)
    {
        _listDistances->setCurrentItem(_listDistances->item(index),QItemSelectionModel::Select);
    }
}

void QMeasureWidget::updateAngles(QStringList list)
{
    bool sel=false;
    int index;
    if (_listAngles->currentItem())
    {
        if (_listAngles->currentItem()->isSelected())
        {
            sel = true;
            index = _listAngles->currentIndex().row();
        }
    }
    _listAngles->clear();
    for (int i=0; i<list.count(); i++)
    {
        _listAngles->addItem(list.at(i));
    }
    if (sel)
    {
        _listAngles->setCurrentItem(_listAngles->item(index),QItemSelectionModel::Select);
    }

}

void QMeasureWidget::updateDihedrals(QStringList list)
{
    bool sel=false;
    int index;
    if (_listDihedrals->currentItem())
    {
        if (_listDihedrals->currentItem()->isSelected())
        {
            sel = true;
            index = _listDihedrals->currentIndex().row();
        }
    }
    _listDihedrals->clear();
    for (int i=0; i<list.count(); i++)
    {
        _listDihedrals->addItem(list.at(i));
    }
    if (sel)
    {
        _listDihedrals->setCurrentItem(_listDihedrals->item(index),QItemSelectionModel::Select);
    }


}

void QMeasureWidget::OnDistanceChange(QListWidgetItem *i)
{
    _listAngles->setCurrentItem(_listAngles->currentItem(),QItemSelectionModel::Deselect);
    _listDihedrals->setCurrentItem(_listDihedrals->currentItem(),QItemSelectionModel::Deselect);
    emit distanceChange(_listDistances->row(i));

}

void QMeasureWidget::OnAngleChange(QListWidgetItem *i)
{
    _listDistances->setCurrentItem(_listDistances->currentItem(),QItemSelectionModel::Deselect);
    _listDihedrals->setCurrentItem(_listDihedrals->currentItem(),QItemSelectionModel::Deselect);
    emit angleChange(_listAngles->row(i));

}

void QMeasureWidget::OnDihedralChange(QListWidgetItem *i)
{
    _listDistances->setCurrentItem(_listDistances->currentItem(),QItemSelectionModel::Deselect);
    _listAngles->setCurrentItem(_listAngles->currentItem(),QItemSelectionModel::Deselect);
    emit dihedralChange(_listDihedrals->row(i));

}

void QMeasureWidget::OnShowDistances()
{
    m_bshowdistances = !m_bshowdistances;
    emit showDistances(m_bshowdistances);
}
