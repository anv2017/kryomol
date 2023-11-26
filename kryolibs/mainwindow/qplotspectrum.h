/*****************************************************************************************
                            qplotspectrum.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QPLOTSPECTRUM_H
#define QPLOTSPECTRUM_H

#include <QWidget>
#include <QHBoxLayout>

#include "fidarray.h"

class QwtPlot;
class QwtPlotCurve;
class QwtPlotZoomer;
class QwtPlotMagnifier;

class QPlotSpectrum: public QWidget
{
    Q_OBJECT

public:
    QPlotSpectrum(QWidget *parent = 0);
    ~QPlotSpectrum();

    enum BaseLinePosition { Top, Middle, Bottom };
    enum SpectrumType {IR, VCD, UV, ECD};
    void PlotSpectrum();
    void InitialPlotSpectrum();
    void SetType(SpectrumType type) {m_type = type;}

public slots:
  void SetData(const std::vector<fidarray>* data,const const fidarray* totaldata,float,float, float);
  void OnIncrease();
  void OnDecrease();
  void OnZoom(bool);
  void ResetZoom();
  void OnPrint();
  void OnBitmapPicture ();
  void OnSVGPicture ();

private:
    QwtPlot* m_plot;
    std::vector<QwtPlotCurve*> m_curves;
    QwtPlotZoomer* m_zoom;
    QwtPlotMagnifier* m_magnifier;
    QHBoxLayout* m_layout;
    const std::vector<fidarray>* m_data;
    const fidarray* m_totaldata;
    float m_max;
    float m_min;
    float m_shift;
    float m_factor;
    SpectrumType m_type;
    BaseLinePosition m_baseline;

};

#endif // QPLOTSPECTRUM_H

