/*****************************************************************************************
                            qnicsgriddialog.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/
#include "qnicsgriddialog.h"
#include "ui_qnicsgriddialog.h"

QNICSGridDialog::QNICSGridDialog(std::pair<kryomol::Coordinate,kryomol::Coordinate>& box, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::QNICSGridDialog)
{
    /*ui->_xFromDoubleSpinBox->setValue(box.first.x);
    ui->_xToDoubleSpinBox->setValue(box.second.x);
    ui->_yFromDoubleSpinBox->setValue(box.first.y);
    ui->_yToDoubleSpinBox->setValue(box.second.y);
    ui->_zFromDoubleSpinBox->setValue(box.first.z);
    ui->_zFromDoubleSpinBox->setValue(box.second.z);
    ui->_zToDoubleSpinBox->setValue(zto);*/
    ui->setupUi(this);
}

QNICSGridDialog::~QNICSGridDialog()
{
    delete ui;
}

std::pair<kryomol::Coordinate,kryomol::Coordinate> QNICSGridDialog::Box() const
{
    std::pair<kryomol::Coordinate,kryomol::Coordinate> box;
    /*box.first.x=ui->_xFromDoubleSpinBox->value();
    box.second.x=ui->_xToDoubleSpinBox->value();
    box.first.y=ui->_yFromDoubleSpinBox->value();
    box.second.y=ui->_yToDoubleSpinBox->value();
    box.first.z=ui->_zFromDoubleSpinBox->value();
    box.second.z=ui->_zToDoubleSpinBox->value();*/
    return box;
}

double QNICSGridDialog::GridResolution() const
{
    return ui->_resolutionSpinBox->value();
}
