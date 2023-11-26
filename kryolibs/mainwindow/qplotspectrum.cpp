/*****************************************************************************************
                            qplotspectrum.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <QVector>
#include <QHBoxLayout>
#include <QMouseEvent>
#include <QPainter>
#include <QFileDialog>
#include <QSvgGenerator>
#include <QPixmap>
#include <QBitmap>
#include <QImage>

#include "qplotspectrum.h"
#include "qwt_plot.h"
#include "qwt_plot_curve.h"
#include "qwt_plot_canvas.h"
#include "qwt_plot_zoomer.h"
#include "qwt_plot_magnifier.h"
#include "qwt_plot_picker.h"
#include "qwt_plot_layout.h"
#include "qwt_scale_div.h"

#include <sstream>
#include <iostream>
#include <vector>

#include "stringtools.h"

using namespace kryomol;

QPlotSpectrum::QPlotSpectrum (QWidget *parent) : QWidget (parent), m_data(nullptr)
{
    m_layout = new QHBoxLayout (this);
    this->setLayout(m_layout);

    m_plot = new QwtPlot(this);
    m_plot->setCanvasBackground(QColor::fromRgb(255,255,255));

    m_layout->addWidget(m_plot);

    m_zoom = new QwtPlotZoomer (m_plot->canvas());
    m_zoom->setEnabled(false);
}

QPlotSpectrum::~QPlotSpectrum()
{}

void QPlotSpectrum::SetData(const std::vector<fidarray>* data, const fidarray* totaldata, float max, float min, float shift)
{
    qDebug() << "number of datasets" << data->size() << Qt::endl;
    for( const auto& x : *data)
    {
        qDebug() << "size of fidarray" << x.size() << Qt::endl;
        for(size_t i=0;i<x.size();++i)
        {
            qDebug() << x[i].real() << Qt::endl;
        }
    }
    m_data = data;
    m_max = max;
    m_min = min;
    m_shift = shift;
    m_totaldata=totaldata;
    PlotSpectrum();
}

void QPlotSpectrum::PlotSpectrum()
{
    if ( m_data == nullptr ) return;

    delete m_zoom;
    delete m_plot;

    m_plot = new QwtPlot(this);
    m_plot->setCanvasBackground(QColor::fromRgb(255,255,255));

    m_layout->addWidget(m_plot);

    m_zoom = new QwtPlotZoomer (m_plot->canvas());
    m_zoom->setEnabled(false);


    const std::vector<fidarray>& datasets=*m_data;

    m_curves.clear();
    m_curves.resize(datasets.size());

    float huestep=0.8/m_curves.size();
    float lowhue=0.1;
    for(auto& c : m_curves)
    {
        c= new QwtPlotCurve();
        c->setStyle(QwtPlotCurve::Lines);
        //Rotate the hue value
        QColor cl;
        cl.setHslF(lowhue,0.4,0.4);
        lowhue+=huestep;
        c->setPen(QPen(cl));
    }



    QVector<double> m_y,m_x;

    double step;
    step = (m_max-m_min)/datasets.front().size();

    for(size_t i=0;i<datasets.size();++i)
    {
        double vx = m_min;
        auto& d=datasets[i];
        auto& c=m_curves[i];
        m_x.clear();
        m_y.clear();
        for (size_t i = 0; i < d.size(); ++i)
        {
            m_x.push_back(vx);
            m_y.push_back(d[i].real());
            vx = vx + step;
        }
        c->setSamples(m_x,m_y);
        c->attach(m_plot);
    }
    m_x.clear();
    m_y.clear();

    //Create total curve
    QwtPlotCurve* ct= new QwtPlotCurve();
    ct->setStyle(QwtPlotCurve::Lines);
    QPen blackpen;
    blackpen.setWidth(blackpen.width()*2);
    blackpen.setColor(Qt::black);
    ct->setPen(blackpen);
    double vx=m_min;
    for(size_t i=0;i<m_totaldata->size();++i)
    {
        m_x.push_back(vx);
        m_y.push_back((*m_totaldata)[i].real());
        vx+=step;
    }
    for(int i=0;i<m_y.size();++i)
    {
        qDebug() << "m_y" << m_y[i] << Qt::endl;
    }
    ct->setSamples(m_x,m_y);
    ct->attach(m_plot);

    QRectF rect = m_curves.front()->boundingRect();

    m_zoom->setZoomBase(rect);
    QString axtitle;
    float baseline;
    switch(m_type)
    {
    case UV:
        axtitle="(eV)";
        baseline=m_zoom->zoomRect().height();
        break;
    case ECD:
        axtitle="(eV)";
        baseline=m_zoom->zoomRect().height()*0.5;
        break;
    case IR:
        axtitle="(cm-1)";
        baseline=m_zoom->zoomRect().height();
        break;
    case VCD:
        axtitle="(cm-1)";
        baseline=m_zoom->zoomRect().height()*0.5;
        break;
    }

    for(auto& c : m_curves)
    {
        c->setBaseline(baseline);
    }
    ct->setBaseline(baseline);

    m_plot->replot();

    update();
}


void QPlotSpectrum::OnIncrease()
{
    const QwtScaleDiv scaleDiv = m_plot->axisScaleDiv(QwtPlot::yLeft);
    switch(m_type)
    {
    case UV:
    case IR:
        m_plot->setAxisScale(QwtPlot::yLeft, scaleDiv.lowerBound(), scaleDiv.lowerBound()+scaleDiv.range()/1.05);
        m_plot->replot();
        break;
    case VCD:
    case ECD:
        m_plot->setAxisScale(QwtPlot::yLeft, scaleDiv.upperBound()-scaleDiv.range()/1.025, scaleDiv.upperBound()+scaleDiv.range()/.025);
        m_plot->replot();
        break;
    }
}


void QPlotSpectrum::OnDecrease()
{
    const QwtScaleDiv scaleDiv = m_plot->axisScaleDiv(QwtPlot::yLeft);
    switch(m_type)
    {
    case UV:
    case IR:
        m_plot->setAxisScale(QwtPlot::yLeft, scaleDiv.lowerBound(), scaleDiv.lowerBound()+scaleDiv.range()*1.05);
        m_plot->replot();
        break;
    case VCD:
    case ECD:
        m_plot->setAxisScale(QwtPlot::yLeft, scaleDiv.upperBound()-scaleDiv.range()*1.025, scaleDiv.upperBound()+scaleDiv.range()*1.025);
        m_plot->replot();
        break;
    }
}


void QPlotSpectrum::OnZoom ( bool b )
{
    if (b)
        m_zoom->setEnabled(true);
    else
    {
        m_zoom->setEnabled(false);
    }
}


void QPlotSpectrum::ResetZoom ()
{
    PlotSpectrum();
}


void QPlotSpectrum::OnPrint()
{
    qDebug() << "reimplement";
    //  QPrinter* pP= new QPrinter ( QPrinter::HighResolution );
    //  pP->setOrientation ( QPrinter::Landscape );

    //  QPrintDialog dialog(pP);

    //  if (dialog.exec())
    //  {
    //    QPainter p ( ( QPrinter* ) pP );
    //    QwtPlotPrintFilter pF;
    //    QRect myrect;
    //    myrect.setRect(0,0,pP->width(),pP->height());
    //    m_plot->print(&p,myrect,pF);
    //    p.end();
    //  }
}


void QPlotSpectrum::OnBitmapPicture ()
{
    qDebug() << "reimplement";
    //    QString filename = QFileDialog::getSaveFileName();

    //    qryomol::StringTokenizer f (filename.toStdString(), ".*");

    //    qryomol::StringTokenizer::iterator it;
    //    int found=-1;
    //    for ( it=f.begin();it!=f.end();++it )
    //    {
    //        qryomol::tolower(*it);
    //        std::string dotext="."+ ( *it);
    //        found=filename.lastIndexOf ( dotext.c_str() );
    //        if ( found != -1 ) break;
    //    }
    //    std::string ext;
    //    if ( found==-1 )
    //    {
    //      filename=filename+".jpeg";
    //      ext="jpeg";
    //    }
    //    else ext= ( *it );

    //    if ( ext=="jpeg" || ext =="jpg" || ext== "png" || ext== "xpm" )
    //    {
    //        QPixmap pix ( QSize(1024,768));

    //        pix.fill(Qt::white);

    //        QPainter p ( &pix );
    //        m_plot->print(&p,QRect(0,0,1024,768));

    //        pix.save(filename,ext.c_str());

    //        p.end();
    //    }

}


void QPlotSpectrum::OnSVGPicture ()
{
    qDebug() << "reimplement";
    //    const QString filename = QFileDialog::getSaveFileName ();
    //    QSvgGenerator picture;
    //    picture.setSize(QSize(800,600));
    //    picture.setFileName(filename);
    //    QPainter p ( &picture );
    //    m_plot->print ( &p, QRect ( 0, 0, 800, 600 ) );
    //    p.end();
}

