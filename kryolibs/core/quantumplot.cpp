/*****************************************************************************************
                            quantumplot.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "quantumplot.h"
#include <qwt_plot_canvas.h>
#include <qwt_plot_layout.h>
#include <qwt_plot_grid.h>
#include <qwt_symbol.h>
//Added by qt3to4:
#include <QMouseEvent>
#include <QToolTip>


QuantumPlot::QuantumPlot ( QWidget* parent ) : QwtPlot ( parent )
{
    setSizePolicy ( QSizePolicy ( QSizePolicy::Preferred,QSizePolicy::Preferred ) );
    connect ( this,SIGNAL ( plotMousePressed ( const QMouseEvent& ) ),this,SLOT ( OnMousePressed ( const QMouseEvent& ) ) );
    //connect(this,SIGNAL(plotMouseMoved(const QMouseEvent& )),this,SLOT(OnMouseMoved(const QMouseEvent& )));

    this->setMouseTracking ( true );

    QwtPlotGrid* grid= new QwtPlotGrid();
    grid->enableX ( true );
    grid->enableY ( true );
    grid->setMajorPen ( QPen ( Qt::gray, 0, Qt::DotLine ) );
    grid->attach ( this );
}


QuantumPlot::~QuantumPlot()
{ }

void QuantumPlot::setupStepMarkers()
{

    const QwtPlotItemList& list=itemList();

    for ( int i=0;i<list.size();++i )
    {

        if ( list.at(i)->rtti() == QwtPlotItem::Rtti_PlotCurve )
        {
            QwtPlotCurve* cv=static_cast<QwtPlotCurve*> ( list.at(i) );
            for ( int j=0;j<cv->dataSize() ;++j )
            {
                QuantumMarker* mrk= new QuantumMarker (this);
                qDebug() << "mrk=" << mrk;
                //         //long marker= insertMarker ( mrk );
                //
                mrk->setSymbol ( new QwtSymbol( QwtSymbol::Diamond, QBrush ( Qt::yellow ), QPen ( Qt::green ), QSize ( 7,7 ) ) );
                mrk->setValue(cv->sample(j));
                mrk->setLineStyle ( QwtPlotMarker::NoLine );
                mrk->setStep ( j );
                mrk->attach ( ( QwtPlot* ) this );
            }
        }
    }

}


void QuantumPlot::mouseMoveEvent(QMouseEvent* e)
{
    const QPoint& pos=e->pos();
    QuantumMarker* mrk= this->closestQuantumMarker ( pos,5 );
    if ( mrk )
    {
        const QwtPlotItemList& list=this->itemList();
        QwtPlotItemList::const_iterator it;
        QString str;
        bool nfirst=false;
        for ( int i=0;i<list.size();++i)
        {

            if ( list.at(i)->rtti() == QwtPlotItem::Rtti_PlotCurve )
            {
                if ( nfirst ) str+="\n";
                str=str+static_cast<QuantumCurve*> ( list.at(i) )->strData ( mrk->step() );
                nfirst=true;
            }
        }
        QRect rect(pos.x(),pos.y(),pos.x()+3,pos.y()+3);

        QToolTip::showText(mapToGlobal(e->pos()),str,this,rect);
    }

}


void QuantumPlot::mousePressEvent (QMouseEvent* e )
{
    QuantumMarker* mrk=closestQuantumMarker ( e->pos(),5 );
    if ( mrk )
    {
        emit selectedPoint ( mrk->step() );
        OnSelectedPoint ( mrk->step() );
    }
    QwtPlot::mousePressEvent ( e );
}

QuantumPlot::QuantumMarker* QuantumPlot::closestQuantumMarker ( const QPoint& pos,int distance )
{
    const QwtPlotItemList& list=itemList();
    QuantumMarker* mrk=NULL;
    QRectF crect=plotLayout()->canvasRect();
    for ( int i=0;i<list.size();++i )
    {
        if ( list.at(i)->rtti() == QwtPlotItem::Rtti_PlotMarker )
        {
            mrk=static_cast<QuantumMarker*> ( list.at(i) );
            int xpos=transform ( QwtPlot::xBottom,mrk->xValue() );
            int ypos=transform ( QwtPlot::yLeft,mrk->yValue() );

            int mxpos=pos.x()-crect.left() ;
            int mypos=pos.y()-crect.top();
            int r=qRound ( sqrt ( ( mxpos-xpos ) * ( mxpos-xpos ) + ( mypos-ypos ) * ( mypos-ypos ) ) );
            if ( r < distance )
            {
                return mrk;
            }
        }
    }
    return NULL;

}

// void QuantumPlot::OnSelectedPoint ( int point )
// {
//   if ( m_point == point ) return;
//   int previous=m_point;
//
//   m_point=point;
//   //change the marker
//   QwtPlotMarkerIterator itm=markerIterator();
//   QwtPlotMarker* m;
//   for ( m=itm.toFirst();itm!=0;m=++itm )
//   {
//     if ( static_cast<QuantumMarker*> ( m )->step() == point )
//     {
//       m->setSymbol ( QwtSymbol ( QwtSymbol::Diamond, red, blue, QSize ( 10,10 ) ) );
//
//     }
//     if ( static_cast<QuantumMarker*> ( m )->step() == previous )
//     {
//       m->setSymbol ( QwtSymbol ( QwtSymbol::Diamond, yellow, green, QSize ( 7,7 ) ) );
//
//     }
//   }
//   replot();
//
// }


void QuantumPlot::OnSelectedPoint ( size_t point )
{
    qDebug() <<"OnSelectedPoint" << point;
    if ( m_point == (int) point ) return;
    int previous=m_point;

    m_point=point;
    //change the marker
    const QwtPlotItemList& list= itemList();

    for ( int i=0;i<(int)list.size();++i )
    {
        if ( list.at(i)->rtti() == QwtPlotItem::Rtti_PlotMarker )
        {
            qDebug() << "reimplement";
            QuantumMarker* m=dynamic_cast<QuantumMarker*> ( list.at(i) );
            if ( m )
            {
                qDebug() << "m=" << m;
                qDebug() << "m->step=" << m->step();
                if ( m->step() == (int) point )
                {
                    qDebug() << "red";
                    //we need to allocate always a new symbol otherwise we have a crush since QwtPlotMarker will destroy the symbol
                    //allocated before
                    m->setSymbol ( new QwtSymbol ( QwtSymbol::Diamond, QBrush ( Qt::red ), QPen ( Qt::blue ), QSize ( 10,10 ) ) );

                }
                if ( m->step() == previous )
                {
                    qDebug() << "yellow";
                    m->setSymbol (  new QwtSymbol( QwtSymbol::Diamond, QBrush ( Qt::yellow ), QPen ( Qt::green ), QSize ( 7,7 ) ) );

                }
            }
        }
    }
    replot();

}
// bool QuantumPlot::setMarkerStep ( long key, int step )
// {
//
//   QwtPlotMarker* m=marker ( key );
//   if ( m )
//   {
//     static_cast<QuantumMarker*> ( m )->setStep ( step );
//     return true;
//   }
//
//   return false;
// }
// long QuantumPlot::closestQuantumMarker ( int x, int y,int& distance )
// {
//   long key=closestMarker ( x,y,distance );
//
//   QuantumMarker* mw= dynamic_cast<QuantumMarker*> ( marker ( key ) );
//   if ( mw ) return key;
//   //Well we dont have a QuantumMarker but a threshold so iterate
//   QwtArray<long> mkeys=markerKeys();
//   QwtArray<long>::Iterator it;
//   QwtPlotMarker* m;
//   double radius=1e5;
//   for ( it=mkeys.begin();it!=mkeys.end();it++ )
//   {
//     if ( dynamic_cast<QuantumMarker*> ( marker ( *it ) ) )
//     {
//       double kx,ky;
//       markerPos ( *it,kx,ky );
//       int r= ( kx-x ) * ( kx-x ) + ( ky-y ) * ( ky-y );
//       if ( radius < r )
//       {
//         radius=r;
//         key=*it;
//       }
//     }
//   }
//
//   distance=sqrt ( radius );
//   return key;
// }


void QuantumPlot::ChangeLimits(int xfirst, int xlast, double ymin,double ymax)
{
    this->setAxisScale ( QwtPlot::xBottom,xfirst,xlast );
    this->setAxisScale( QwtPlot::yLeft,ymin,ymax);
    this->replot();
}
