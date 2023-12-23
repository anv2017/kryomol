/*****************************************************************************************
                            frequency.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef FREQUENCY_H
#define FREQUENCY_H

#include "coordinate.h"
#include <vector>

struct Frequency
{
    Frequency() : x{0},y{0},z{0},w{0} { }
    Frequency(float fx, float fy,float fz=0.0,float fw=0.0) { x=fx; y=fy; z=fz; w=fw; }
    float x,y,z,w;
};

class Spectralline
{
public:
  Spectralline() { x=y0=y1=y2=y3=0.0f; m_solventshift=0; }
  Spectralline(float fx, float fy0, float fy1, float fy2, float fy3) : x(fx), y0(fy0), y1(fy1), y2(fy2), y3(fy3) {}

  float x; float y0; float y1; float y2; float y3;
  kryomol::Coordinate& ElectricDipole() { return m_electricdipole; }
  const kryomol::Coordinate& ElectricDipole() const { return m_electricdipole; }
  kryomol::Coordinate& VelocityDipole() { return m_velocitydipole; }
  const kryomol::Coordinate& VelocityDipole() const { return m_velocitydipole; }
  kryomol::Coordinate& MagneticDipole() { return m_magneticdipole; }
  const kryomol::Coordinate& MagneticDipole() const { return m_magneticdipole; }
  float RotatoryStrengthLength() const { return m_rotatorystrenghtlength; }
  void SetRotatoryStrengthLength(float r) { m_rotatorystrenghtlength=r; }
  float RotatoryStrengthVelocity() const { return m_rotatorystrenghtvelocity; }
  void SetRotatoryStrengthVelocity(float r) { m_rotatorystrenghtvelocity=r; }
  void SetSolventShift(float shift) { m_solventshift=shift; }
  float SolventShift() const { return m_solventshift; }
  std::vector<float> Coefficient() const { return m_coefficient; }
  void SetCoefficient(std::vector<float>  r) { m_coefficient=r; }
  std::vector<int> OrbitalI() const { return m_orbitali; }
  void SetOrbitalI(std::vector<int> i) { m_orbitali=i; }
  std::vector<int> OrbitalJ() const { return m_orbitalj; }
  void SetOrbitalJ(std::vector<int> j) { m_orbitalj=j; }


private:
  kryomol::Coordinate m_electricdipole;
  kryomol::Coordinate m_magneticdipole;
  kryomol::Coordinate m_velocitydipole;
  float m_rotatorystrenghtvelocity;
  float m_rotatorystrenghtlength;
  float m_solventshift;
  std::vector<float> m_coefficient;
  std::vector<int> m_orbitali;
  std::vector<int> m_orbitalj;

};

#endif
