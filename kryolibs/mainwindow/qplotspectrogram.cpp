/*****************************************************************************************
                            qplotspectrogram.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#if QT_VERSION >= 0x040000
#include <qprintdialog.h>
#endif
#include <qwt_color_map.h>
#include <qwt_plot_spectrogram.h>
#include <qwt_scale_widget.h>
#include <qwt_scale_draw.h>
#include <qwt_plot_zoomer.h>
#include <qwt_plot_panner.h>
#include <qwt_plot_layout.h>

#include <iostream>
#include <QDebug>

#include "qplotspectrogram.h"
#include "density.h"
#include "mathtools.h"

using namespace kryomol;

class Zoom: public QwtPlotZoomer
{
public:
    Zoom(QWidget *canvas):
        QwtPlotZoomer(canvas)
    {
        setTrackerMode(AlwaysOn);
    }

    virtual QwtText trackerText(const QPointF &pos) const
    {
        QColor bg(Qt::white);
#if QT_VERSION >= 0x040300
        bg.setAlpha(200);
#endif

        QwtText text = QwtPlotZoomer::trackerTextF(pos);
        text.setBackgroundBrush( QBrush( bg ));
        return text;
    }

};

class SpectrogramData: public QwtRasterData
{
public:
    SpectrogramData(int nx, int ny, float dx, float dy, float max, float min, const D2Array<float>& planedata) : QwtRasterData()
        //QwtRasterData(QRectF(-((nx-1)*dx)/2, -((ny-1)*dy)/2, (nx-1)*dx, (ny-1)*dy)), m_planedata(planedata)
    {
        //qDebug() << "SpectrogramData: " << -((nx-1)*dx)/2 << " " << -((ny-1)*dy)/2 << " " << (nx-1)*dx << " " << (ny-1)*dy << " " << m_planedata.NRows() << " " << m_planedata.NColumns() << endl;
        //qDebug() << "Rect: " << this->boundingRect().height() << " " << this->boundingRect().width() << " " << this->boundingRect().topLeft() << " " << this->boundingRect().topRight() << " " << this->boundingRect().bottomLeft() << " " << this->boundingRect().bottomRight() << endl;
        m_nx = nx;
        m_ny = ny;
        m_dx = dx;
        m_dy = dy;
        m_max = max;
        m_min = min;
    }

    virtual QwtRasterData *copy() const
    {
        return new SpectrogramData(m_nx, m_ny, m_dx, m_dy, m_max, m_min, m_planedata);
    }

    virtual QwtInterval interval() const
    {
        return QwtInterval(m_min, m_max);
    }

    virtual double value(double x, double y) const
    {
        size_t i = (size_t) floor(m_nx/2 + x/m_dx + 0.5);
        size_t j = (size_t) floor(m_ny/2 + y/m_dy + 0.5);

        return m_planedata(i,j);
    }

private:
    int m_nx;
    int m_ny;
    float m_dx;
    float m_dy;
    float m_max;
    float m_min;
    D2Array<float> m_planedata;
};

QPlotSpectrogram::QPlotSpectrogram(kryomol::Density *density, QWidget *parent): QwtPlot(parent), m_density(density)
{
    m_spectrogram = new QwtPlotSpectrogram();
    m_spectrogram->attach(this);
    m_zoomer = new Zoom(canvas());
    m_spectrogram->attach(this);

    m_plane = QPlotSpectrogram::XY;
    m_level = 100000;
}

void QPlotSpectrogram::FillSpectrogram()
{
    if (!m_density->DensityMatrix().Empty())
    {
        m_zoomer->deleteLater();        

        if (m_level==100000)
            SetAxisValue(0.0);


        qDebug() << "SPECTROGRAM: " << endl;
        qDebug() << "m_density: " << m_density->Nx() << " " << m_density->Ny() << " " << m_density->Nz() << " " << m_density->Dx() << " " << m_density->Dy() << " " << m_density->Dz() << " m_level: " << m_level << endl;

        QwtScaleWidget *xAxis = axisWidget(QwtPlot::xBottom);;
        QwtScaleWidget *yAxis = axisWidget(QwtPlot::yLeft);
        D2Array<float> planedata;
        switch (m_plane)
        {
            case QPlotSpectrogram::XY:
                planedata = m_density->DensityMatrix().XY(m_level);
                qDebug() << "Planedata: " << planedata.NRows() << " " << planedata.NColumns() << endl;
                m_spectrogram->setData( new SpectrogramData(m_density->Nx(),m_density->Ny(),m_density->Dx(),m_density->Dy(),m_density->DensityMatrix().Max(),m_density->DensityMatrix().Min(),planedata));
                xAxis->setTitle("X");
                setAxisScale(QwtPlot::xBottom, -((m_density->Nx()-1)/2)*m_density->Dx(), ((m_density->Nx()-1)/2)*m_density->Dx());
                yAxis->setTitle("Y");
                setAxisScale(QwtPlot::yLeft, -((m_density->Ny()-1)/2)*m_density->Dy(), ((m_density->Ny()-1)/2)*m_density->Dy());
            break;
            case QPlotSpectrogram::XZ:
                planedata = m_density->DensityMatrix().XZ(m_level);
                m_spectrogram->setData(new SpectrogramData(m_density->Nx(),m_density->Nz(),m_density->Dx(),m_density->Dz(),m_density->DensityMatrix().Max(),m_density->DensityMatrix().Min(),planedata));
                xAxis->setTitle("X");
                setAxisScale(QwtPlot::xBottom, -((m_density->Nx()-1)/2)*m_density->Dx(), ((m_density->Nx()-1)/2)*m_density->Dx());
                yAxis->setTitle("Z");
                setAxisScale(QwtPlot::yLeft, -((m_density->Nz()-1)/2)*m_density->Dz(), ((m_density->Nz()-1)/2)*m_density->Dz());
            break;
            case QPlotSpectrogram::YZ:
                planedata = m_density->DensityMatrix().YZ(m_level);
                m_spectrogram->setData(new SpectrogramData(m_density->Ny(),m_density->Nz(),m_density->Dy(),m_density->Dz(),m_density->DensityMatrix().Max(),m_density->DensityMatrix().Min(),planedata));
                xAxis->setTitle("Y");
                setAxisScale(QwtPlot::xBottom, -((m_density->Ny()-1)/2)*m_density->Dy(), ((m_density->Ny()-1)/2)*m_density->Dy());
                yAxis->setTitle("Z");
                setAxisScale(QwtPlot::yLeft, -((m_density->Nz()-1)/2)*m_density->Dz(), ((m_density->Nz()-1)/2)*m_density->Dz());

            break;
        }
        //qDebug() << "Spectrogram Data: " << m_spectrogram->data()->boundingRect().bottomLeft() << " " << m_spectrogram->data().boundingRect().bottomRight() << " " << m_spectrogram->data().boundingRect().topLeft() << " " << m_spectrogram->data().boundingRect().topRight() << endl;


        //qDebug() << "Maximum: " << m_density->DensityMatrix().Max() << " Minimum: " << m_density->DensityMatrix().Min() << endl;
        QList<double> contourLevels;
        for ( double level = m_density->DensityMatrix().Min() + (m_density->DensityMatrix().Max()-m_density->DensityMatrix().Min())/41; level < m_density->DensityMatrix().Max(); level += (m_density->DensityMatrix().Max()-m_density->DensityMatrix().Min())/41)
            contourLevels.push_back(level);
        m_spectrogram->setContourLevels(contourLevels);

        qDebug() << " check ";
        QwtInterval interval=m_spectrogram->data()->interval(Qt::XAxis);
        double range= interval.maxValue()-interval.minValue();
        QwtLinearColorMap colorMap(QColor(143,0,255), QColor(255,0,0));
        colorMap.addColorStop((0.0-interval.minValue())/(4*range), QColor(143, 0, 255));
        colorMap.addColorStop(2*(0.0-interval.minValue())/(4*range), QColor(0, 0, 255));
        colorMap.addColorStop(3*(0.0-interval.minValue())/(4*range), QColor(0, 255, 0));
        colorMap.addColorStop((0.0-interval.minValue())/range, QColor(0,0,0));
        colorMap.addColorStop(5*(0.0-interval.minValue())/(4*range), QColor(255, 255, 0));
        colorMap.addColorStop(6*(0.0-interval.minValue())/(4*range), QColor(255, 127, 0));

        m_spectrogram->setColorMap(&colorMap);

        // A color bar on the right axis
        QwtScaleWidget *rightAxis = axisWidget(QwtPlot::yRight);
        rightAxis->setTitle("Intensity");
        rightAxis->setColorBarEnabled(true);
        QwtInterval yrinterval=m_spectrogram->data()->interval(Qt::YAxis);
        rightAxis->setColorMap(yrinterval, const_cast<QwtColorMap*>(m_spectrogram->colorMap())) ;
        setAxisScale(QwtPlot::yRight, yrinterval.minValue(), yrinterval.maxValue() );

        enableAxis(QwtPlot::xBottom);
        enableAxis(QwtPlot::yLeft);
        enableAxis(QwtPlot::yRight);

        qDebug() << "1----" << endl;

        plotLayout()->setAlignCanvasToScales(true);
        replot();

        qDebug() << "2----" << endl;

        // LeftButton for the zooming
        // MidButton for the panning
        // RightButton: zoom out by 1
        // Ctrl+RighButton: zoom out to full size

        m_zoomer = new Zoom(canvas());
    #if QT_VERSION < 0x040000
        m_zoomer->setMousePattern(QwtEventPattern::MouseSelect2,
            Qt::RightButton, Qt::ControlButton);
    #else
        m_zoomer->setMousePattern(QwtEventPattern::MouseSelect2,
            Qt::RightButton, Qt::ControlModifier);
    #endif
        m_zoomer->setMousePattern(QwtEventPattern::MouseSelect3,
            Qt::RightButton);

        qDebug() << "3----" << endl;

        QwtPlotPanner *panner = new QwtPlotPanner(canvas());
        panner->setAxisEnabled(QwtPlot::yRight, false);
        panner->setMouseButton(Qt::MidButton);

        qDebug() << "4----" << endl;

        // Avoid jumping when labels with more/less digits
        // appear/disappear when scrolling vertically

        const QFontMetrics fm(axisWidget(QwtPlot::yLeft)->font());
        QwtScaleDraw *sd = axisScaleDraw(QwtPlot::yLeft);
        sd->setMinimumExtent( fm.width("100.00") );

        qDebug() << "5----" << endl;

        const QColor c(Qt::darkBlue);
        m_zoomer->setRubberBandPen(c);
        m_zoomer->setTrackerPen(c);

        this->setAutoReplot(true);

        qDebug() << "6----" << endl;
    }
}

void QPlotSpectrogram::showContour(bool b)
{
    m_spectrogram->setDisplayMode(QwtPlotSpectrogram::ContourMode, b);
    FillSpectrogram();
}

void QPlotSpectrogram::showSpectrogram(bool b)
{
    m_spectrogram->setDisplayMode(QwtPlotSpectrogram::ImageMode, b);
    m_spectrogram->setDefaultContourPen(b ? QPen() : QPen(Qt::NoPen));
    FillSpectrogram();
}

void QPlotSpectrogram::printSpectrogram()
{
    qDebug() << "reimplement";
//    QPrinter printer;
//    printer.setOrientation(QPrinter::Landscape);
//#if QT_VERSION < 0x040000
//    printer.setColorMode(QPrinter::Color);
//#if 0
//    printer.setOutputFileName("/tmp/spectrogram.ps");
//#endif
//    if (printer.setup())
//#else
//#if 0
//    printer.setOutputFileName("/tmp/spectrogram.pdf");
//#endif
//    QPrintDialog dialog(&printer);
//    if ( dialog.exec() )
//#endif
//    {
//        print(printer);
//    }
}

void QPlotSpectrogram::SetAxisValue(double value)
{
    switch (m_plane)
    {
        case QPlotSpectrogram::XY:
            m_level = floor((m_density->Nz())/2+value/m_density->Dz());
            break;
        case QPlotSpectrogram::XZ:
            m_level = floor((m_density->Ny())/2+value/m_density->Dy());
            break;
        case QPlotSpectrogram::YZ:
            m_level = floor((m_density->Nx())/2+value/m_density->Dx());
            break;
    }
}
