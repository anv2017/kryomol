/*****************************************************************************************
                            irspectrum.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef IRSPECTRUM_H
#define IRSPECTRUM_H
#include <vector>

#include "sinusoid.h"
#include "frequency.h"
#include "qplotspectrum.h"

class IRSpectrum
{
public:
  IRSpectrum();
  ~IRSpectrum();
  void SetFrequencies(const std::vector< std::vector<Frequency> >& v,double scale=1.0);
  void CalculateSpectrum();
  bool WriteJCampDX();
  void CopyData();
  void SetFile(const char* file);
  const std::vector<fidarray>&   GetData() const { return m_data; }
  const fidarray& GetTotalData() const { return m_totaldata; }
  const fidarray& GetExperimentalData() const { return m_expdata; }
  float Max() { return m_max; }
  float Min() { return m_min; }
  float Shift() { return m_shift;}
  int NPoints() { return m_npoints;}
  float LineWidth() const { return m_linewidth; }
  void SetLimits(float& max, float& min) { m_max=max; m_min=min;}// CalculateSpectrum(); }
  void SetAutoLimits();
  void SetLineWidth(float lw);
  void SetShift(float shift);
  void SetNPoints (int n) {m_npoints = n;}
  void SetType(QPlotSpectrum::SpectrumType type);
  QPlotSpectrum::SpectrumType GetType() const { return m_spectrumtype; }
private:
  float GetIntensityAt(float lambda,const std::vector<Sinusoid>& sv);
protected:
  std::string m_file;
  /** a vector of frequencies for each frame */
  std::vector< std::vector<Frequency> > m_frequencysets;
private:
  /** a vector of sinusoids for each frame*/
  std::vector< std::vector<Sinusoid> > m_sinusoidsets;
  std::vector<fidarray> m_data;
  fidarray m_expdata;
  fidarray m_totaldata;
  float m_max;
  float m_min;
  std::string m_title;
  QPlotSpectrum::SpectrumType m_spectrumtype;
  /** Linewidth in cm-1*/
  float m_linewidth;
  float m_shift;
  int m_npoints;
  /** conformational weights*/
  std::vector<float> m_weights;
};

#endif
