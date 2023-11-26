/*****************************************************************************************
                            qlogtable.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "qlogtable.h"
#include "qtablewidget.h"
#include "qheaderview.h"

#include "mathtools.h"

#include <math.h>

QLogTable::QLogTable(QWidget* parent)
    : QTableWidget(parent)
{
    m_order=Qt::AscendingOrder;
}


QLogTable::~QLogTable()
{}


void QLogTable::sortItems ( int col)
{
    int i, j, k;
    QString temp;

    for (i = (rowCount() - 1); i >= 0; i--)
    {
      for (j = 1; j <= i; j++)
      {
          if(m_order == Qt::AscendingOrder)
        {
          if (item(j-1,col)->text().toDouble() > item(j,col)->text().toDouble())
          {
            for (k = 0; k < columnCount(); k++)
            {
                temp = item(j-1,k)->text();
                item(j-1,k)->setText(item(j,k)->text());
                item(j,k)->setText(temp);
            }
            temp = verticalHeaderItem(j-1)->text();
            verticalHeaderItem(j-1)->setText(verticalHeaderItem(j)->text());
            verticalHeaderItem(j)->setText(temp);
          }
        }
        else
        {
          if (item(j-1,col)->text().toDouble() < item(j,col)->text().toDouble())
          {
              for (k = 0; k < columnCount(); k++)
              {
                  temp = item(j-1,k)->text();
                  item(j-1,k)->setText(item(j,k)->text());
                  item(j,k)->setText(temp);
              }
              temp = verticalHeaderItem(j-1)->text();
              verticalHeaderItem(j-1)->setText(verticalHeaderItem(j)->text());
              verticalHeaderItem(j)->setText(temp);
          }

        }

      }
    }

    if (m_order == Qt::AscendingOrder)
    {
        m_order = Qt::DescendingOrder;
        horizontalHeader()->setSortIndicator(col,Qt::DescendingOrder);
        horizontalHeader()->setSortIndicatorShown(true);
    }
    else
    {
        m_order = Qt::AscendingOrder;
        horizontalHeader()->setSortIndicator(col,Qt::AscendingOrder);
        horizontalHeader()->setSortIndicatorShown(true);
    }

}

