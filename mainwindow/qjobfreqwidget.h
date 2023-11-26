/*****************************************************************************************
                            qjobfreqwidget.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QJOBFREQWIDGET_H
#define QJOBFREQWIDGET_H

#include <QWidget>

#include "qjobwidget.h"
#include "qfreqwidget.h"
#include "qplotspectrum.h"

namespace kryomol
{
class World;

};


class QJobFreqWidget : public QJobWidget
{
    Q_OBJECT

public:
    QJobFreqWidget(const QString& file, kryomol::World* world, QWidget *parent = 0);
    ~QJobFreqWidget();

private:
    void Init();

private slots:
    void OnShowSpectrum (bool bshow);
    void OnIRTypeChanged(QPlotSpectrum::SpectrumType);

private:
    QString m_file;
    QFreqWidget* m_freqwidget;


};

#endif // QJOBFREQWIDGET_H
