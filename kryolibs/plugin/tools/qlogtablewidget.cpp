/*****************************************************************************************
                            qlogtablewidget.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <math.h>

#include <QTextStream>
#include <QMenu>
#include <QMouseEvent>
#include <QClipboard>
#include <QApplication>
#include <QVariant>

#include "mathtools.h"
#include "qlogtablewidget.h"

using namespace kryomol;

template <class X> class QTableWidgetFloatingNumItem : public QTableWidgetItem
{
public:
    QTableWidgetFloatingNumItem(X data, char format = 'g', int precision = 6 ) : _data(data), _format(format), _precision(precision)
    {

    }
    QVariant data(int role ) const
    {
        if ( role == Qt::DisplayRole )
        {

            QString str;
            str.setNum(_data,_format,_precision);
            return str;
        }
        else
        {
            //setData(role,);
            return QTableWidgetItem::data(role);
        }

    }
    bool operator < (const QTableWidgetItem& other) const
    {
        const QTableWidgetFloatingNumItem* nOther=static_cast<const QTableWidgetFloatingNumItem*>(&other);
        return ( this->_data < nOther->_data );
    }
private:
    X _data;
    char _format;
    int _precision;

};

class QTableWidgetIntItem : public QTableWidgetItem
{
public:
    QTableWidgetIntItem(int data ) : _data(data)
    {

    }
    QVariant data(int role ) const
    {
        if ( role == Qt::DisplayRole )
        {
            return QString::number(_data);
        }
        else
        {
            //setData(role,);
            return QTableWidgetItem::data(role);
        }

    }
    bool operator < (const QTableWidgetItem& other) const
    {
        const QTableWidgetIntItem* nOther=static_cast<const QTableWidgetIntItem*>(&other);
        return ( this->_data < nOther->_data );
    }
private:
    int _data;

};

class QLogTableWidgetPrivate
{
  public:
    QLogTableWidgetPrivate() {}
    ~QLogTableWidgetPrivate() {}
  public:
    QMenu* m_menu;
    Qt::SortOrder m_order;
    bool m_fabssorting;
    bool m_copyvheader;
    bool m_copyhheader;
};

QLogTableWidget::QLogTableWidget ( QWidget* parent ) : QTableWidget ( parent )
{
  Init();
}

QLogTableWidget::QLogTableWidget ( int rows, int columns,QWidget* parent ) : QTableWidget ( rows,columns,parent )
{
  Init();
}

//void QLogTableWidget::setAbsoluteValueSorting ( bool b )
//{
//  _d->m_fabssorting=b;
//}

void QLogTableWidget::setCopyHorizontalHeader ( bool b )
{
  _d->m_copyhheader=b;
}

void QLogTableWidget::setCopyVerticalHeader ( bool b )
{
  _d->m_copyvheader=b;
}
void QLogTableWidget::Init()
{
  _d=new QLogTableWidgetPrivate();
  _d->m_menu= new QMenu ( this );
  QAction* act=_d->m_menu->addAction ( "Copy Contents to clipboard" );
  connect ( act,SIGNAL ( triggered() ),this,SLOT ( OnCopyContents() ) );
  _d->m_order=Qt::DescendingOrder;
  _d->m_fabssorting=false;
  _d->m_copyvheader=false;
  _d->m_copyhheader=false;
}
QLogTableWidget::~QLogTableWidget()
{
  delete _d;
}

//void QLogTableWidget::sortColumn ( int column, Qt::SortOrder order )
//{

//  qDebug() << "column pressed is " << column << endl;
//  bubbleSort ( column,_d->m_order == Qt::AscendingOrder );

//  if ( _d->m_order == Qt::AscendingOrder )
//  {
//    _d->m_order=Qt::DescendingOrder;
//  }
//  else _d->m_order=Qt::AscendingOrder;
//}

//void QLogTableWidget::bubbleSort ( int col, bool ascending )
//{
//  int i, j;

//  int ncols=columnCount() +1; //include header
//  std::vector< D1Array<double> > values ( rowCount() );
//  for ( i=0;i<rowCount();++i )
//    values[i].Initialize ( ncols );
//  std::vector< D1Array<double>* > pv ( rowCount() );
//  D1Array<double>* pkv;
//  for ( i=0;i< rowCount();++i )
//    for ( j=0;j< ncols; ++j )
//    {
//      if ( j == 0 )
//      {
//        values[i] ( j ) =verticalHeaderItem ( i )->text().toDouble();
//      }
//      else
//      {
//        values[i] ( j ) =item ( i,j-1 )->text().toDouble();
//      }
//      pv[i]=&values[i];
//    }

//  for ( i = ( rowCount() - 1 ); i >= 0; i-- )
//  {
//    for ( j = 1; j <= i; j++ )
//    {
//      if ( ascending )
//      {
//        if ( _d->m_fabssorting )
//        {
//          if ( fabs ( ( ( *pv[j-1] ) ) ( col+1 ) ) < fabs ( ( *pv[j] ) ( col+1 ) ) )
//          {
//            pkv=pv[j];
//            pv[j]=pv[j-1];
//            pv[j-1]=pkv;

//          }
//        }
//        else
//        {
//          if ( ( ( *pv[j-1] ) ) ( col+1 ) < ( *pv[j] ) ( col+1 ) )
//          {
//            pkv=pv[j];
//            pv[j]=pv[j-1];
//            pv[j-1]=pkv;

//          }
//        }
//      }
//      else
//      {
//        if ( _d->m_fabssorting )
//        {
//          if ( fabs ( ( *pv[j-1] ) ( col+1 ) ) > fabs ( ( *pv[j] ) ( col+1 ) ) )
//          {
//            pkv=pv[j];
//            pv[j]=pv[j-1];
//            pv[j-1]=pkv;
//          }
//        }
//        else
//        {
//          if ( ( *pv[j-1] ) ( col+1 ) > ( *pv[j] ) ( col+1 ) )
//          {
//            pkv=pv[j];
//            pv[j]=pv[j-1];
//            pv[j-1]=pkv;
//          }
//        }

//      }

//    }
//  }

//  QString s;
//  for ( i=0;i<rowCount();++i )
//  {
//    for ( j=0;j<ncols;++j )
//    {
//      if ( j == 0 )
//      {
//        s.sprintf ( "%d", static_cast<int> ( ( *pv[i] ) ( j ) ) );
//        verticalHeaderItem ( i )->setText ( s );
//      }
//      else
//      {
//        s.sprintf ( "%.3f", ( *pv[i] ) ( j ) );
//        item ( i,j-1 )->setText ( s );
//      }
//    }
//  }
//}


void QLogTableWidget::mousePressEvent ( QMouseEvent* e )
{
  if ( e->button() == Qt::RightButton )
  {
    _d->m_menu->exec ( QCursor::pos() );
  }
  QTableWidget::mousePressEvent ( e );
}

void QLogTableWidget::OnCopyContents()
{
  QString str;
  QTextStream stream ( &str );

  QTableWidgetItem* it;

   if (_d->m_copyhheader )
   {
     for(int j=0;j<columnCount();++j)
     {
       it=this->horizontalHeaderItem(j);
       stream << it->text() << '\t';
     }
   }
   stream << endl;

  for ( int i=0;i<rowCount();++i )
  {

    if ( _d->m_copyvheader) 
    {
      it=this->verticalHeaderItem(i);
      stream << it->text() << '\t';
    }
    
    for ( int j=0;j<columnCount();++j )
    {

      it=item ( i,j );
      stream << it->data(Qt::DisplayRole).toString() << '\t';
    }
    stream << endl;
  }
  QClipboard* clipboard = QApplication::clipboard();
  if ( clipboard )
  {
    clipboard->setText ( str );
  }
}


void QLogTableWidget::setNum(int row, int col, float num, char format,int precision)
{
    QTableWidgetFloatingNumItem<float>* item= new QTableWidgetFloatingNumItem<float>(num,format,precision);
    item->setFlags ( Qt::ItemIsEnabled | Qt::ItemIsSelectable );
    this->setItem(row,col,item);
}

void QLogTableWidget::setNum(int row, int col, double num, char format,int precision)
{
    QTableWidgetFloatingNumItem<double>* item= new QTableWidgetFloatingNumItem<double>(num,format,precision);
    item->setFlags ( Qt::ItemIsEnabled | Qt::ItemIsSelectable );
    this->setItem(row,col,item);
}

void QLogTableWidget::setNum(int row, int col, int num)
{
    QTableWidgetIntItem* item = new QTableWidgetIntItem(num);
    item->setFlags ( Qt::ItemIsEnabled | Qt::ItemIsSelectable );
    this->setItem(row,col,item);

}


