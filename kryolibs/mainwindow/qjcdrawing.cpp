/*****************************************************************************************
                            qjcdrawing.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <iostream>

#include <qpainter.h>
#include <qpixmap.h>
#include <qpaintdevice.h>
#include <qmessagebox.h>
#include <qcursor.h>
#include <qmenu.h>

#include <qwt_scale_draw.h>
#include <qwt_scale_engine.h>

#include <QRubberBand>
#include "qjcdrawing.h"
//Added by qt3to4:
#include <QPaintEvent>
#include <QResizeEvent>
#include <QMouseEvent>
#include "utilities.h"
//#include "tokenizer.h"
#include "cursors.h"
#include "stringtools.h"

const float increasepass=1.3f;
const int scalecoarse=10;
const int scalefine=10;
QJCDrawing::QJCDrawing ( QWidget *parent)
    : QWidget ( parent)
{
    setAccessibleName ( "JCDrawing" );
    m_scale=  new QwtScaleDraw();
    m_scalepalette= new QPalette();
    m_scalepalette->setColor ( QPalette::Text, QColor ( 0,0,0 ) );
    m_scalepalette->setColor ( QPalette::Foreground, QColor ( 0,0,0 ) );
    m_scaleengine = new QwtLinearScaleEngine();

    setBackgroundRole ( QPalette::NoRole );
    m_color.setRgb ( 0,0,0 );
    m_bkcolor.setRgb ( 255,255,255 );
    m_factor=1.0f;
    //no data in construnction
    m_data=NULL;

    //tracker for zoom
    m_tracker=new QRubberBand ( QRubberBand::Rectangle,this );

    //build up menu
    m_optionsmenu =new QMenu ( "spectrum options menu", this);
    m_optionsmenu->addAction ( "Set Zoom Limits as Spectrum Limits",this,SLOT ( OnSetLimits() ) );
    m_optionsmenu->addAction ( "Copy SVG to Clipboard",this,SLOT ( OnCopy() ) );
    m_orientation=RightLeft;
    m_baselineposition=Bottom;
}


QJCDrawing::~QJCDrawing()
{}



void QJCDrawing::resizeEvent ( QResizeEvent* e )
{
    UpdateDrawing();
}

void QJCDrawing::paintEvent ( QPaintEvent* e )
{
    QPixmap pix ( this->rect().size() );
    QPainter p ( &pix );
    Draw ( &p,this->rect() );
    p.end();
    p.begin ( this );
    p.drawPixmap ( this->rect().topLeft(), pix );

}

void QJCDrawing::Draw ( QPainter* pDC, QRect rect, bool bprint )
{
    m_rect=rect;
    pDC->fillRect ( rect,m_bkcolor );


    if ( m_data )
    {
        if ( bprint )
            DrawPrintData ( pDC );
        else
            DrawData ( pDC );
        switch ( m_orientation )
        {
        case RightLeft:
            m_scale->setScaleDiv ( m_scaleengine->divideScale ( m_min+ ( m_max-m_min ) * ( 1-m_zoom.fleft ),m_min+ ( m_max-m_min ) * ( 1-m_zoom.fright ) ,scalecoarse,scalefine,0.0 ) );

            m_scale->draw ( pDC,*m_scalepalette);//->active() );
            break;
        case LeftRight:
            m_scale->setScaleDiv ( m_scaleengine->divideScale ( m_min+ ( m_max-m_min ) * ( m_zoom.fleft ),m_min+ ( m_max-m_min ) * ( m_zoom.fleft ),scalecoarse,scalefine ) );
            m_scale->draw ( pDC, *m_scalepalette );

            break;
        }
    }

}

void QJCDrawing::UpdateDrawing ( bool autoscale /*=false*/ )
{
    UpdateDrawing ( autoscale,this->rect() );
}
void QJCDrawing::UpdateDrawing ( bool autoscale, const QRect& rect )
{

    //we need to resample the data set to the widht of the drawing rect
    //const QRect& rect=this->rect();
    const fidarray& data=*m_data;
    float max=-1e12;
    float min=1e12;
    int l=m_zoom.firstpoint;
    int r=m_zoom.lastpoint;

    if ( autoscale )
    {
        for ( int i=l;i<r;i++ )
        {
            if ( data[i].real() > max ) max=data[i].real();
            if ( data[i].real() < min ) min=data[i].real();
        }
        float fbaseline;
        switch ( m_baselineposition )
        {
        case Middle:
            fbaseline=0.5;
            break;
        default:
            fbaseline=0.8;
            break;

        }
        int baseline=qRound ( rect.top() +rect.height() *fbaseline );
        m_factor=- ( rect.top() +5-baseline ) /max;
        max=-1e12;
        min=1e12;
    }


    m_points.clear();
    float xscalefactor=static_cast<float> ( rect.width() ) /static_cast<float> ( r-l );
    float value=0;
    int x=-1;
    float fbaseline;
    switch ( m_baselineposition )
    {
    case Middle:
        fbaseline=0.5;
        break;
    default:
        fbaseline=0.8;
        break;

    }
    int baseline=qRound ( rect.top() +rect.height() *fbaseline );

    for ( int i=l;i<r;++i )
    {
        if ( data[i].real() > max ) max=data[i].real();
        if ( data[i].real() < min ) min=data[i].real();
        int newx=qRound ( ( i-l ) *xscalefactor );
        if ( newx!= x )
        {
            m_points.push_back ( QPoint ( x,baseline-qRound ( m_factor*min ) ) );
            m_points.push_back ( QPoint ( x,baseline-qRound ( m_factor*max ) ) );
            value=0;
            x=newx;
            max=-1e12;
            min=1e12;
        }

    }

    m_scale->move( QPoint(rect.left(),rect.top()+rect.height()*0.9));
    m_scale->setLength(rect.width());

}

void QJCDrawing::DrawData ( QPainter* p )
{
    if ( !m_points.empty() )
        p->drawPolyline(&m_points[0],m_points.size());
}

void QJCDrawing::DrawPrintData ( QPainter* p )
{
#warning maybe not necessary
}

void QJCDrawing::SetData ( fidarray* data,float max, float min, float shift, int npoints )
{

    m_data=data;
    m_min=min;
    m_max=max;
    m_shift = shift;
    m_npoints = npoints;
    m_zoom.SetPoints ( data->size() );
    UpdateDrawing ( true );
    update();


}


void QJCDrawing::OnIncrease()
{
    m_factor*=increasepass;
    UpdateDrawing();
    update();
}


void QJCDrawing::OnDecrease()
{
    m_factor/=increasepass;
    UpdateDrawing();
    update();
}

void QJCDrawing::OnPrint()
{
    qDebug() << "reimplement";
    //  QPrinter* pP= new QPrinter ( QPrinter::HighResolution );
    //  pP->setOrientation ( QPrinter::Landscape );

    //  QPrintDialog dialog(pP);

    //  if (dialog.exec()) //( pP->setup() )
    //  {
    //    QPainter p ( ( QPrinter* ) pP );
    //    QRect myrect;
    //    myrect.setRect(0,0,pP->width(),pP->height());

    //    UpdateDrawing ( true,myrect );
    //    Draw ( &p,myrect );
    //    UpdateDrawing(true);
    //    update();
    //  }


}

void QJCDrawing::Export ( const QString& s,const QString& filter )
{

    QMessageBox kk;
    kk.setText ( filter );
    kk.exec();
    kryomol::StringTokenizer f (filter.toStdString(), ".*");

    kryomol::StringTokenizer::iterator it;
    int found=-1;
    for ( it=f.begin();it!=f.end();++it )
    {
        kryomol::tolower(*it);
        std::string dotext="."+ ( *it);
        found=s.lastIndexOf ( dotext.c_str() );
        if ( found != -1 ) break;
    }
    QString ms;
    std::string ext;
    if ( found==-1 )
    {
        ms=s+"."+QString ( f.front().c_str() );
        ext=f.front();
    }
    else ext= ( *it );

    if ( ext=="svg" )
    {
        GetPicture().save ( ms,ext.c_str() );
        return;
    }
    std::cerr << ext << std::endl;
    if ( ext=="jpeg" || ext =="jpg" || ext== "png" || ext== "xpm" )
    {
        kryomol::toupper ( ext );
        GetImage().save ( ms,ext.c_str() );
        return;
    }
}

QPicture QJCDrawing::GetPicture()
{
    QPicture pic;
    QPainter p;
    p.begin ( &pic );
    QRect rect ( 0,0,1024,768 );
    UpdateDrawing ( false,rect );
    Draw ( &p,rect,true );
    p.end();

    UpdateDrawing();
    update();
    return pic;
}

QImage QJCDrawing::GetImage()
{
    QRect rect ( 0,0,1024,768 );
    QPixmap pix ( rect.size() );
    QPainter p ( &pix );
    UpdateDrawing ( false,rect );
    Draw ( &p,rect );
    UpdateDrawing();
    update();
    return pix.toImage();

}


void QJCDrawing::mousePressEvent ( QMouseEvent* e )
{
    switch ( e->button() )
    {
    case Qt::LeftButton:

        switch ( m_visualmode )
        {
        case eZOOM:
            m_tracker->setGeometry(QRect(e->pos(),QSize()));
            m_tracker->show();
            break;
            /*   case eINFO:
               m_grabbed=true;
               linetracker->begin(e->pos());
               break;*/

        default:
            break;
        }

        break;
    case Qt::RightButton:
        m_optionsmenu->exec ( QCursor::pos() );
        break;
    default:
        break;
    }

}

void QJCDrawing::mouseMoveEvent ( QMouseEvent* e )
{
    switch ( m_visualmode )
    {
    case eZOOM:
    {
        QPoint topleft=m_tracker->geometry().topLeft();
        m_tracker->setGeometry(QRect(topleft,e->pos()).normalized());
    }
        break;
    default:
        break;
    }
}

void QJCDrawing::mouseReleaseEvent ( QMouseEvent* e )
{

    if ( e->button() ==Qt::LeftButton )
    {
        QRect Rect,rct;
        float ratio,offsetl,offsetr;
        switch ( m_visualmode )
        {
        case eZOOM:
        {
            QRect rct=m_tracker->geometry();
            m_tracker->hide();
            if ( rct.width() >6 )
            {
                Rect=this->rect();
                ratio= ( float ) rct.width() / ( float ) Rect.width();
                if ( rct.left() <= Rect.left() ) offsetl=0;
                else
                    offsetl=fabs ( ( ( float ) ( rct.left()-Rect.left() ) ) / ( float ) Rect.width() );
                if ( rct.right() >= Rect.right() ) offsetr=1;
                else
                    offsetr=fabs ( ( ( float ) ( rct.right()-Rect.left() ) ) / ( float ) Rect.width() );

                SetZoom ( offsetl, offsetr );
                UpdateDrawing();
                update();
            }
        }
        default:
            break;
        }
    }
}

/*!
    \fn QJCDrawing::SetZoom(float l, float r)
 */
void QJCDrawing::SetZoom ( float r, float o )
{
    int left;
    int right;

    switch ( m_orientation )
    {
    case RightLeft:
        left=m_zoom.firstpoint;
        right=m_zoom.lastpoint;
        m_zoom.firstpoint=qRound ( left+r* ( right -left ) );
        m_zoom.lastpoint=qRound ( left+o* ( right-left ) );
        break;
    case LeftRight:
        left=m_zoom.lastpoint;
        right=m_zoom.firstpoint;
        m_zoom.lastpoint=qRound ( left-r* ( left-right ) );
        m_zoom.firstpoint=qRound ( left-o* ( left-right ) );
        break;
    }

    float nleft=r* ( m_zoom.fright-m_zoom.fleft );
    float nright=o* ( m_zoom.fright-m_zoom.fleft );
    m_zoom.fright=nright+m_zoom.fleft;
    m_zoom.fleft+=nleft;

}

void QJCDrawing::SetCursor()
{

    switch ( m_visualmode )
    {
    case eZOOM:
        this->setCursor ( QCursor ( QPixmap ( cursors::zoom ),14,14 ) );
        break;
    default:
        this->setCursor ( QCursor ( Qt::ArrowCursor ) );
        break;
    }
}


void QJCDrawing::ResetZoom()
{
    m_zoom.reset ( m_data->size() );
    UpdateDrawing ( true );
    update();
}

void QJCDrawing::OnSetLimits()
{
    float left;
    float right;
    switch ( m_orientation )
    {
    case RightLeft:
        left=m_min+ ( m_max-m_min ) * ( 1-m_zoom.fleft );
        right=m_min+ ( m_max-m_min ) * ( 1-m_zoom.fright );
        emit limits ( left,right );
        break;
    case LeftRight:
        left=m_min+ ( m_max-m_min ) * ( m_zoom.fleft );
        right=m_min+ ( m_max-m_min ) * ( m_zoom.fright );
        emit limits ( right,left );
        break;
    }
    ResetZoom();

}


void QJCDrawing::OnCopy()
{
#warning reimplement
    /*
  QPicMimeSource* s= new QPicMimeSource ( &GetPicture() );
  s->copy();
*/
}

void QJCDrawing::OnSnapshot()
{
#warning reimplement
    /*
  QDir dir=QDir::home();
  QString sfile;

  KFileDialog* fd= new KFileDialog (
                       QString::null,
                       "*.svg|SVG vectorialformat (*.svg)\n*.png *.PNG|PNG raster format (*.png)\n*.jpeg *.jpg *.JPEG *.JPG|JPEG raster format (*.jpeg *.jpg)",
                       this,"",true );


  if ( fd->exec() == QDialog::Accepted )
    sfile = fd->selectedFile();


  if ( sfile==QString::null )
    return;
  QString filter=fd->currentFilter();

  Export ( sfile,filter );
*/
}


void QJCDrawing::OnZoom ( bool b )
{

    if ( b ) SetVisualMode ( eZOOM );
    else SetVisualMode ( eNORMAL );
}

void QJCDrawing::SetBaseLinePosition ( BaseLinePosition b )
{
    m_baselineposition=b;
    UpdateDrawing ( true );
    update();
}
