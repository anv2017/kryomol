/*****************************************************************************************
                            coordinate.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "coordinate.h"
#include <cmath>


using namespace kryomol;

Coordinate::Coordinate() : m_xyz(3,0.)
{}

Coordinate::Coordinate(float cx,float cy,float cz) : m_xyz(3)
{
  x()=cx;
  y()=cy;
  z()=cz;
}

Coordinate::~Coordinate()
{}

//vectorial product
Coordinate Coordinate::operator ^(const Coordinate& c) const
{

  Coordinate v;
  v.x()=this->y()*c.z()-this->z()*c.y();
  v.y()=this->z()*c.x()-this->x()*c.z();
  v.z()=this->x()*c.y()-this->y()*c.x();
  return v;
}


float Coordinate::operator *(const Coordinate& b) const
{
   return this->x()*b.x()+this->y()*b.y()+this->z()*b.z();
}

void Coordinate::operator-=(const Coordinate& b)
{
  this->x()-=b.x();
  this->y()-=b.y();
  this->z()-=b.z();
}

void Coordinate::operator+=(float f)
{
  this->x()+=f;
  this->y()+=f;
  this->z()+=f;
}

void Coordinate::operator+=(const Coordinate& b)
{
  this->x()+=b.x();
  this->y()+=b.y();
  this->z()+=b.z();
}


void Coordinate::operator*=(float f)
{
  this->x()*=f;
  this->y()*=f;
  this->z()*=f;
}

void Coordinate::operator/=(float f)
{
  this->x()/=f;
  this->y()/=f;
  this->z()/=f;
}

Coordinate Coordinate::operator*(float f) const
{
  Coordinate a=(*this);
  a*=f;
  return a;
}

Coordinate Coordinate::operator/(float f) const
{
  Coordinate a=(*this);
  a/=f;
  return a;
}

Coordinate Coordinate::operator-(const Coordinate& b) const
{
  Coordinate a=*this;
  a-=b;
  return a;

}

Coordinate Coordinate::operator+(const Coordinate& b) const
{
  Coordinate a=*this;
  a+=b;
  return a;

}

float Coordinate::Distance(const Coordinate& a, const Coordinate& b)
{
  float x1=a.x()-b.x();
  float y1=a.y()-b.y();
  float z1=a.z()-b.z();
  return sqrt( x1*x1+y1*y1+z1*z1 );
}

float Coordinate::Angle(const Coordinate& a, const Coordinate& b)
{
  float proj=a*b/(a.Norm()*b.Norm()) ;
  if ( proj < -1 ) proj=-1;
  if (proj > 1) proj=1;
  return acos(proj);
}

float Coordinate::Norm() const
{
  return sqrt(x()*x()+y()*y()+z()*z());
}

void Coordinate::Normalize()
{
#ifdef __GNUC__
#warning performance
#endif
  float norm=this->Norm();
  x()/=norm;
  y()/=norm;
  z()/=norm;
}

float Coordinate::ScalarProduct(const Coordinate& a, const Coordinate& b)
{
    return a.x()*b.x()+a.y()*b.y()+a.z()*b.z();
}

float Coordinate::Dihedral(const Coordinate& a, const Coordinate& b, const Coordinate& c, const Coordinate& d)
{
  Coordinate vec1,vec2,vec3,v1,v2;


  vec1=a-b;
  vec2=d-c;
  vec3=c-b;

  /*vec3.x=1.0;
  vec3.y=0.0;
  vec3.z=0.0;
  vec2.x=0.0;
  vec2.y=1.0;
  vec2.z=0.0;*/
  //v1=vec3**vec2
  v1=vec3^vec2;

  //v2=vec1*vec3
  v2=vec3^vec1;

  float ang=Angle(v1,v2);


  //OK lets determine the sign
  Coordinate dp=vec1^vec2;
  float proj=dp*vec3;
  if(proj<0)
    ang*=-1.0f;
  return ang;
}

float Coordinate::GetAngle(const Coordinate& a1, const Coordinate& a2, const Coordinate& a3,bool degrees /*=false*/)
{

  float angle=Coordinate::Angle(a1-a2,a3-a2);

  if(degrees)
    return angle*180/M_PI;
  else
    return angle;
}

float Coordinate::GetDihedral(const Coordinate& a1, const Coordinate& a2, const Coordinate& a3, const Coordinate& a4,bool degrees /*=false*/)
{

  float d= Coordinate::Dihedral(a1,a2,a3,a4);

  if(degrees)
    return d*180/M_PI;
  else
    return d;
}

Coordinate Coordinate::MiddlePoint(const Coordinate& a, const Coordinate& b)
{
  Coordinate c;
  c.x()=0.5*(a.x()+b.x());
  c.y()=0.5*(a.y()+b.y());
  c.z()=0.5*(a.z()+b.z());
  return c;
}

//Rotate using Quaternions
Coordinate Coordinate::RotAroundAxis(const Coordinate& c, const Coordinate& axisorigin, const Coordinate& axisend, float theta)
{
  Coordinate axis=axisend-axisorigin;
  std::cout << "axis is " << axis << std::endl;
  std::cout << "theta is " << theta << std::endl;
  float kk=(sin(theta/2)/axis.Norm());
  Coordinate v1=axis*kk;

  Quaternion q(cos(theta/2),v1);
  Quaternion q1=q.Invert();

  Quaternion q2(0,c-axisorigin);

  Quaternion q3=q^q2^q1;
  std::cout << "q3 " << q3 << std::endl;

  return axisorigin+q3.V();


}

