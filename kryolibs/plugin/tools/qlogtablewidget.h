/*****************************************************************************************
                            qlogtablewidget.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QLOGTABLEWIDGET_H
#define QLOGTABLEWIDGET_H

#include <QTableWidget>

#include "export.h"

class QLogTableWidgetPrivate;

/** @brief The base class for qryomol plugins
    A QTableWidget derived class for output of floating point values.
    Allows to copy values to spreadsheets and fast sorting by value*/


class KRYOMOL_API QLogTableWidget : public QTableWidget
{
    Q_OBJECT

public:
    QLogTableWidget(QWidget* parent=0);
    QLogTableWidget(int rows, int columns,QWidget* parent=0);
    ~QLogTableWidget();
    ///** order the columns according to absolute value (false by default*/
    //   void setAbsoluteValueSorting(bool b);
    /** copy vertical header along with contents*/
    void setCopyVerticalHeader(bool b);
    /** copy vertical header along with contents*/
    void setCopyHorizontalHeader(bool b);
    /** set the cell i,j to number n)*/
    void setNum(int row, int col, float num, char format='g',int precision=6);
    void setNum(int row, int col, double num, char format='g',int precision=6);
    void setNum(int row, int col, int num);

protected:
    void mousePressEvent(QMouseEvent* e);

private:
    //void bubbleSort(int column,bool ascending);
    void Init();

public slots:
    //void sortColumn(int column,Qt::SortOrder order= Qt::AscendingOrder);

private slots:
    void OnCopyContents();

private:
    QLogTableWidgetPrivate* _d;
};

#endif
