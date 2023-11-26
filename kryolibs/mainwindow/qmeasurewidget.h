/*****************************************************************************************
                            qmeasurewidget.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QMEASUREWIDGET_H
#define QMEASUREWIDGET_H

#include <QWidget>

#include "ui_qmeasurewidgetbase.h"


class QMeasureWidget : public QWidget, private Ui::QMeasureWidgetBase
{
    Q_OBJECT

public:
    explicit QMeasureWidget(QWidget *parent = 0);
    ~QMeasureWidget();

    void Init();
    void updateDistances(QStringList list);
    void updateAngles(QStringList list);
    void updateDihedrals(QStringList list);

public slots:
    void OnClearAll();
    void OnWriteDistance(QString&);
    void OnWriteAngle(QString&);
    void OnWriteDihedral(QString&);

signals:
    void clearAll();
    void distanceChange(int);
    void angleChange(int);
    void dihedralChange(int);
    void showDistances(bool);

private slots:
    void OnDistanceChange(QListWidgetItem*);
    void OnAngleChange(QListWidgetItem*);
    void OnDihedralChange(QListWidgetItem*);
    void OnShowDistances();

private:
    bool m_bshowdistances;

};

#endif // QMEASUREWIDGET_H
