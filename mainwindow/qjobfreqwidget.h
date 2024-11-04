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
#include "qplotspectrum.h"

namespace kryomol
{
class World;

};

class QIRWidget;
class QFreqWidget;

class QJobFreqWidget : public QJobWidget
{
    Q_OBJECT

public:
    QJobFreqWidget(const QString& file, QWidget *parent = 0);
    ~QJobFreqWidget();

    void InitWidgets();

private slots:
    void OnShowSpectrum (bool bshow);
    void OnIRTypeChanged(QPlotSpectrum::SpectrumType);

private:
    QString m_file;
    QFreqWidget* m_freqwidget;
    QIRWidget* m_irwidget;


};

#endif // QJOBFREQWIDGET_H
