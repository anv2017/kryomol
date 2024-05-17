/*****************************************************************************************
                            uvspectrum.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef UVSPECTRUM_H
#define UVSPECTRUM_H

/**
A class for simulation of UV and ECD spectra
*/

#include <vector>

#include "kryomolcore_export.h"
#include "frequency.h"
#include "fidarray.h"
#include "sinusoid.h"
#include "qplotspectrum.h"

class KRYOMOLCORE_EXPORT UVSpectrum
{
public:
  enum formalism { length=0, velocity };
  UVSpectrum();
  ~UVSpectrum();
  void SetLines(std::vector< std::vector<Spectralline> >& v,double scale=1.0);
  void CalculateSpectrum();
  bool WriteJCampDX();
  void CopyData();
  void SetFile(const char* file);
  const std::vector<fidarray>*   GetData() const { return &m_data; }
  const fidarray* GetTotalData() const { return &m_totaldata; }
  const fidarray& GetExperimentalData() const { return m_expdata; }
  float Max() { return m_max; }
  float Min() { return m_min; }
  float Shift() { return m_shift;}
  int NPoints() { return m_npoints;}
  void SetAutoLimits();
  void SetLimits(float max, float min) { m_max = max; m_min = min;}
  void SetLineWidth(float lw);
  void SetShift (float sh);
  void SubstractSolventShift(bool b);
  void SetNPoints (int n) {m_npoints = n;}
  float LineWidth() const { return m_linewidth; }
  void SetType(QPlotSpectrum::SpectrumType type);
  QPlotSpectrum::SpectrumType GetType() const { return m_spectrumtype; }
  formalism Formalism() const { return m_formalism; }
  void SetFormalism(formalism f);
  void SetEnantiomer(bool b);
  bool Enantiomer() const { return m_benantiomer; }
  /** set use boltzmann weighthing for averaged curve*/
  void SetBoltzmannWeighting(bool b) { m_boltzw=b; }
  /** true if using boltzmann weighthing for averaged curve*/
  bool BoltzmannWeighting() const { return m_boltzw; }
  void SetPopulations(const std::vector<double>* p) { m_populations=p; }



private:
  float GetIntensityAt(float lambda,const std::vector<Sinusoid>& sdd);
  void RecalculateX();
protected:
  std::string m_file;
  std::vector< std::vector<Spectralline> > m_linesets;
private:
  std::vector< std::vector<Sinusoid> > m_sinusoidsets;
  /** there will be a curve for each of the n conformations*/
  std::vector<fidarray> m_data;
  /** sum of the n curves*/
  fidarray m_totaldata;
  /** experimental curve*/
  fidarray m_expdata;
  float m_max;
  float m_min;
  std::string m_title;
  QPlotSpectrum::SpectrumType m_spectrumtype;
  float m_linewidth;
  float m_shift;
  int m_npoints;
  formalism m_formalism;
  bool m_benantiomer;
  bool m_bsubstractsolventshift;
  std::vector<float> m_weights;
  bool m_boltzw;
  const std::vector<double>* m_populations;
};

#endif
