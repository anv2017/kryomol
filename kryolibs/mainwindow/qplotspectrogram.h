/*****************************************************************************************
                            qplotspectrogram.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QPLOTSPECTROGRAM_H
#define QPLOTSPECTROGRAM_H

#include <qwt_plot.h>
#include <qwt_plot_spectrogram.h>
#include <qwt_plot_zoomer.h>

namespace kryomol
{
  class Density;
}

class QPlotSpectrogram: public QwtPlot
{
    Q_OBJECT

public:
    enum Plane{XY,XZ,YZ};
    QPlotSpectrogram(kryomol::Density *density, QWidget * = NULL);
    void FillSpectrogram();
    void SetAxisPlane(Plane plane) {m_plane=plane;}
    void SetAxisValue(double value);

public slots:
    void showContour(bool on);
    void showSpectrogram(bool on);
    void printSpectrogram();

private:
    QwtPlotSpectrogram *m_spectrogram;
    QwtPlotZoomer* m_zoomer;
    kryomol::Density *m_density;
    Plane m_plane;
    size_t m_level;
};

#endif // QPLOTSPECTROGRAM_H
