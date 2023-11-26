/*****************************************************************************************
                            qfreqwidget.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QFREQWIDGET_H
#define QFREQWIDGET_H

#include <vector>

#include <QWidget>

#include "ui_qfreqwidgetbase.h"
#include "frequency.h"
#include "irspectrum.h"
#include "qplotspectrum.h"
#include "world.h"

class IRSpectrum;
class fidarray;
class QLogTableWidget;

class QFreqWidget : public QWidget,  public IRSpectrum, private Ui::QFreqWidgetBase
{
  Q_OBJECT
public:
    QFreqWidget(kryomol::World* world,const QString& file, QWidget* parent=0);

    ~QFreqWidget();

    /** Init frequency table for conformation with index fidx*/
    void InitTable(size_t fidx);
    /** get a copy of the frequencies for each conformation and pass it to IRSpectrum*/
    void InitFrequencies();
private:
    QFreqWidget(QWidget* parent=0,const char* name=0);
public slots:
    void OnShowSpectrum();
    void OnSetLimits(float,float);
    void OnDisableControls(bool);
signals:
    void showspectrum(bool);
    void data(const std::vector<fidarray>& data,float,float,float);
    void mode(int,int);
    void movie(bool);
    void distort(int);
    void reset();
    void Type(QPlotSpectrum::SpectrumType);
protected:
  QSize sizeHint() const { return QSize(350,400); }
private slots:
  void OnWriteJCampDX();
  void OnCopy();
  void OnChangeLimits();
  void OnResetLimits();
  void OnScaling();
  void OnAnimate(int);
  void OnSetLineWidth(double );
  void OnSetShift(int);
  void OnSetNPoints(int);
  void OnDistort(int);
  void OnResetDistortions();
  void OnSpectrumTypeChanged(int);
  void OnTableSelection(int );
private:
    bool m_bshowspectrum;
    std::vector< std::vector<Frequency> > m_frequencysets;
    std::vector<int> m_distortframes;
    int m_activemode;
    int m_npoints;
    QLogTableWidget* m_freqtable;
    kryomol::World* m_world;

};

#endif
