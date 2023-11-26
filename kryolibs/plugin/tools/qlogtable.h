/*****************************************************************************************
                            qlogtable.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QLOGTABLE_H
#define QLOGTABLE_H

#include <qtableview.h>
#include <qtablewidget.h>
/**
A table with sorting
*/
class QLogTable : public QTableWidget
{
public:
    QLogTable(QWidget* parent=0);

    ~QLogTable();

    void sortItems(int);

private:
  Qt::SortOrder m_order;

};

#endif
