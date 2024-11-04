/*****************************************************************************************
                            qjobspwidget.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QJOBSPWIDGET_H
#define QJOBSPQWIDGET_H

#include <QWidget>

#include "qjobwidget.h"
#include "qfreqwidget.h"
#include "qplotspectrum.h"

namespace kryomol
{
class World;

};


class QJobSpWidget : public QJobWidget
{
    Q_OBJECT

public:
    QJobSpWidget(QWidget *parent = 0);
    ~QJobSpWidget();
    void InitWidgets();

private:
    void Init();

private:
    QString m_file;

};

#endif // QJOBSPWIDGET_H
