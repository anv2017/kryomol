/*****************************************************************************************
                            qconvwidget.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QCONVWIDGET_H
#define QCONVWIDGET_H

#include <iostream>

#include <QWidget>
#include <QSlider>
#include <QCheckBox>

#include "kryomolcore_export.h"
#include "threshold.h"

/**
This class contais QwtPlot widget for info on convergence of calculations
*/

class QuantumPlot;
class Convergence
{
public:
  Convergence(size_t size);
  ~Convergence();
  double* Energies() { return m_energies; }
  double* S2() { return m_S2; }
  double* RMSForces()  {return m_rmsforces; }
  double* MaximumForces() { return m_maximumforces; }
  double* RMSDisplacements() { return m_rmsdisplacements; }
  double* MaximumDisplacements() { return m_maximumdisplacements; }
  double MaxEnergy();
  double MinEnergy();

  size_t size() const { return m_size; }
  friend std::ostream& operator << (std::ostream& s, Convergence& conv);
private:
  double* m_energies;
  double* m_S2;
  double* m_rmsforces;
  double* m_maximumforces;
  double* m_rmsdisplacements;
  double* m_maximumdisplacements;
  size_t m_size;
private:
  double GetMaxArray(double* array);
  double GetMinArray(double* array);

};

class QryoMol;
class QwtPlot;
class QDoubleEditBox;
class QSlider;
class QPushButton;

class KRYOMOLCORE_EXPORT QConvWidget : public QWidget
{
  Q_OBJECT
public:
    QConvWidget(QWidget *parent = 0, bool gradients=true, bool displacements=true );

  ~QConvWidget();
  void SetNData(size_t size);
  double* GetEnergies() { return m_convdata->Energies(); }
  double* GetS2() { return m_convdata->S2(); }
  double* GetRMSForces() { return m_convdata->RMSForces(); }
  double* GetMaximumForces() { return m_convdata->MaximumForces(); }
  double* GetRMSDisplacements() { return m_convdata->RMSDisplacements(); }
  double* GetMaximumDisplacements() { return m_convdata->MaximumDisplacements(); }
  int GetForceScale() { return m_scaleslider->value();}

  void SetupCurves();
  void SetThreshold(const Threshold& thr) { m_threshold=thr; }
  void SetEnergyLevel(const QString& level) {  m_energylevel=level; }
signals:
  void selectedPoint(size_t );
  void forcescale(float );
  void showforces(bool);
public slots:
  void OnSelectedPoint(size_t );
  void OnShowForces();
private slots:
  void OnChangeLimits(int, int);
  void OnForceSliderChanged(int );
private:
  QuantumPlot* m_energy;
  QuantumPlot* m_gradient;
  QuantumPlot* m_displacement;
  double* m_xaxis;
  double* m_yenergies;
  double* m_yrmsforces;
  double* m_ymaximumforces;
  double* m_yrmsdisplacements;
  double* m_ymaximumdisplacements;
  size_t m_ndata;
  QString m_energylevel;
  Convergence* m_convdata;
  Threshold m_threshold;
  QDoubleEditBox* m_limitsbox;
  QSlider* m_scaleslider;
  QCheckBox* m_showforce;
  bool m_bshowforces;
protected:
  QSize sizeHint() const { return QSize(350,250); }

};

#endif
