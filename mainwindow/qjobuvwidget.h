/*****************************************************************************************
                            qjobuvwidget.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QJOBUVWIDGET_H
#define QJOBUVWIDGET_H

#include <QObject>

#include "qjobwidget.h"
#include "quvwidget.h"
#include "qplotspectrum.h"

namespace kryomol
{
class World;

};

class QIRWidget;
class QJobUVWidget : public QJobWidget
{

    Q_OBJECT

public:
    QJobUVWidget(const QString& file, QWidget *parent = 0);
    ~QJobUVWidget();
    void InitWidgets();

private:
    void Init();

private slots:
    void OnShowUVSpectrum (bool bshow);
    void OnUVTypeChanged(QPlotSpectrum::SpectrumType);
    void OnFrameChanged(size_t fidx);

private:
    bool m_beta;
    QString m_file;
    QUVWidget* m_uvwidget;
    QIRWidget* m_irwidget;

};

#endif // QJOBUVWIDGET_H
