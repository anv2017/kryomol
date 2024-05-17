/*****************************************************************************************
                            quvwidget.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QUVWIDGET_H
#define QUVWIDGET_H

#include <QWidget>
#include <vector>
#include "ui_quvwidgetbase.h"
#include "uvspectrum.h"
#include "kryomolcore_export.h"
#include "qplotspectrum.h"
#include "transitionchange.h"

/**
Widget containing controls for diplay of UV and ECD spectra
*/


class QLogTableWidget;
class QVBoxLayout;

namespace kryomol {
class World;
}

class KRYOMOLCORE_EXPORT QUVWidget : public QWidget, public UVSpectrum, private Ui::QUVWidgetBase
{
  Q_OBJECT


public:
    QUVWidget(const kryomol::World* world, const QString& file, QWidget* parent=0);
    ~QUVWidget();
    void InitTable(size_t fidx);
    void SetCheckableTransitionChanges();
    void SetCheckableTransitionCoefficients();    
    void SetTransitionChanges(const std::vector< std::vector< std::vector<kryomol::TransitionChange> > >& transitions)
    { m_transitions=transitions;}
    void SetBeta(bool b) {m_beta = b;}
    void UpdatePlot();
signals:
  void showspectrum(bool);
  void data(const std::vector<fidarray>* d,const fidarray* td,float,float,float);
  void type(QPlotSpectrum::SpectrumType);
  void showelectricdipole(bool );
  void showmagneticdipole(bool );
  void showvelocitydipole(bool );
  void showtransitionchanges(bool);
  void transitionselected(int );
  void showtransition(int );
  void showdensities(int );
  void offshowtransitions(bool);
  void setactivetransition(int);
private slots:
  //void OnSortItems(int);
  void OnWriteJCampDX();
  void OnCopy();
  void OnShowSpectrum();
  void OnSpectrumTypeChanged(int);
  void OnSetLineWidth(double);
  void OnSetShift(int);
  void OnSetNPoints(int);
  void OnChangeLimits();
  void OnResetLimits();
  void OnSetLimits(float, float);
  void OnTransitionSelected(int );
  void OffShowTransitionChanges(bool);
  void OnShowTransitionChanges(bool);
  void OnShowTransitionCoefficients(bool);
  void OnShowDensityChanges(bool);
  void OnShowElectricDipole(bool);
  void OnShowMagneticDipole(bool);
  void OnShowVelocityDipole(bool);
  void OnCalculateEnantiomer(bool );
  void OnFormalismChanged(int );
  void OnSolventShift(bool b);
  void OnBoltzmannCheckBox(bool b);
private:
  void InitTransitionTable(int );

  //void OnRenderOptionsChange();

private:
  QUVWidget(QWidget* parent=0, const char* name=0);
private:
  bool m_bshowspectrum;
  bool m_bshowcoefficients;
  bool m_bshowtransitions;
  bool m_bshowdensities;
  bool m_bshowelectricdipole;
  bool m_bshowmagneticdipole;
  bool m_bshowvelocitydipole;
  bool m_bsolventshift;
  bool m_beta;
  /** a vector holding coefficients for orbital transitions for each spectral line for each conformation*/
  std::vector< std::vector< std::vector<kryomol::TransitionChange> > > m_transitions;
  QLogTableWidget* m_uvTable;
  QVBoxLayout* m_tablelay;
  const kryomol::World* m_world;


};

#endif
