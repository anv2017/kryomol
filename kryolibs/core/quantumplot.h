/*****************************************************************************************
                            quantumplot.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QUANTUMPLOT_H
#define QUANTUMPLOT_H

#include <qwt_plot.h>
#include <qwt_plot_marker.h>
#include <qwt_plot_curve.h>
#include <QMouseEvent>

/**
A derived QwtPlot class which plots markers on all optimization points and emits a signal with the number of the step clicked
*/


class QuantumPlot : public QwtPlot
{
    Q_OBJECT

  public:
    QuantumPlot ( QWidget* parent );

    ~QuantumPlot();

    void setupStepMarkers();
    void ChangeLimits(int first, int last, double ymin,double ymax);
  signals:
    void selectedPoint ( size_t );
  public slots:
    void OnSelectedPoint ( size_t point );
  protected:
    void mousePressEvent ( QMouseEvent* e );
    void mouseMoveEvent( QMouseEvent* e );

  class QuantumMarker : public QwtPlotMarker
    {
      public:
        QuantumMarker ( QwtPlot* parent ) : QwtPlotMarker() {}
        void setStep ( int step ) { m_step=step; }
        int step() const { return m_step; }
      private:
        int m_step;
    };


    QSize sizeHint() const { return QSize ( 250,200 ); }

    int m_point;
public:
    QuantumMarker* closestQuantumMarker ( const QPoint& pos, int distance );
};

class QuantumCurve : public QwtPlotCurve
{
  public:
    QuantumCurve ( QwtPlot* parent, const QString& title=QString::null )
        : QwtPlotCurve ( title ) , m_label(QString::null) {}
    void setExtraData ( double* z ) { m_extradata=z; }
    double z ( int i ) { return m_extradata[i]; }
    void setExtraLabel ( const QString& label ) { m_label=label; }
    const QString& extraLabel() const { return m_label; }
    virtual QString strData ( int i )
    {
      QString str;
      str.sprintf ( "%.6f",sample(i).y() );
      str=extraLabel() +"=" + str;
      return str;
    };
  private:
    double* m_extradata;
    QString m_label;

};

class EnergyCurve :  public QuantumCurve
{
  public:
    EnergyCurve ( QwtPlot* parent, const QString& title=QString::null )
        : QuantumCurve ( parent,title ) {}
    QString strData ( int i )
    {
      QString str;
      str.sprintf ( "%.5f\nS2=%.3f",sample( i ).y(),z ( i ) );
      str=extraLabel()  + "=" +str;
      return str;
    }
};



#endif
