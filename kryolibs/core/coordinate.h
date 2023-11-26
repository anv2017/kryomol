/*****************************************************************************************
                            coordinate.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef COORDINATE_H
#define COORDINATE_H

#include <iostream>
#include <cmath>
#include "mathtools.h"
#include "coreexport.h"

namespace kryomol
{

/** @brief representaion of a 3D point

This class manage representation and basic geometric
manipulation of 3D vectors. Data are stored as a single precision IMArray<float> vector*/
class KRYOMOLCORE_API Coordinate
{
public:
  /** Buid (0,0,0) vector*/
  Coordinate();
  /** Build (cx, cy,cz ) vector*/
  Coordinate(float cx,float cy,float cz);
  /** vectorial product*/
  Coordinate operator^(const Coordinate& b) const;
  /** scalar (dot) product*/
  float operator*(const Coordinate& b) const;
  /** difference of vectors*/
  void operator-=(const Coordinate& b);
  /** sum of vector and scalar*/
  void operator +=(float f);
  /** sum of two vectors*/
  void operator +=(const Coordinate& b);
  /**  multiplication with a scalar*/
  void operator *=(float f);
  /** division by a scalar number*/
  void operator /=(float f);
  /** vector scalar product*/
  Coordinate operator *(float f) const;
  /** vector scalar division*/
  Coordinate operator /(float f) const;
  /** vector-vector difference*/
  Coordinate operator -(const Coordinate& b) const;
  /** vector vector addition*/
  Coordinate operator+(const Coordinate& b) const;
  /** @return x coordinate*/
  float& x() { return m_xyz(0); }
    /** @return x coordinate*/
  const float& x() const { return m_xyz(0); }
    /** @return y coordinate*/
  float& y() { return m_xyz(1); }
 /** @return y coordinate*/
  const float& y() const { return m_xyz(1); }
 /** @return z coordinate*/
  float& z() { return m_xyz(2); }
 /** @return z coordinate*/
  const float& z() const { return m_xyz(2); }
  /** return the @see IMArray representation of the coordinate*/ 
  operator const D1Array<float>& () const { return m_xyz; }
  /** return the @see IMArray representation of the coordinate*/ 
  operator D1Array<float>& () { return m_xyz; }
  /** return the scalar product of two coordinate vectors*/
  static float ScalarProduct(const Coordinate& a, const Coordinate& b);
  /** */
  ~Coordinate();
  /** @return distance between 3D points a and b*/
  static float Distance(const Coordinate& a, const Coordinate& b);
  /** @return the norm of this vector*/
  float Norm() const;
  /** normalize the vector*/
  void Normalize();
  /** @return angle between vectors a and b*/
  static  float Angle(const Coordinate& a, const Coordinate& b);
  /** @return dihedral angle between vectors a,b,c,d*/
  static float Dihedral(const Coordinate& a, const Coordinate& b, const Coordinate& c, const Coordinate& d);
  static float GetAngle(const Coordinate& a, const Coordinate& b, const Coordinate& c, bool degrees);
  /** @return dihedral angle between vectors a,b,c,d*/
  static float GetDihedral(const Coordinate& a, const Coordinate& b, const Coordinate& c, const Coordinate& d, bool degrees);
  /** @return (b-a)/2 */
  static  Coordinate MiddlePoint(const Coordinate& a, const Coordinate& b);
  
  /** Rotate coordinate c theta radians clockwise around axis */
  static Coordinate RotAroundAxis(const Coordinate& c, const Coordinate& axisorigin, const Coordinate& axisend,float angle);

private:
  D1Array<float> m_xyz;
};

/** @brief quaternion representation
  
  an utlilty class for management of quaternions
  each quaternion is represented as combination of a Coordinate v and an scalar w
*/
  class KRYOMOLCORE_API Quaternion
  {
  public:
    /** Built quaternion(0,Coordinate(0,0,0))*/
    Quaternion() { _w=_v.x()=_v.y()=_v.z()=0; }
    /** Built quaternion (scalar,Coordiante(x,y,z))*/
    Quaternion(double scalar,const Coordinate& c) { _w=scalar; _v=c; }  
    /** @return 3d coordinate part v of the quaternion*/
    Coordinate& V() { return _v; }
    /** @return 3d coordinate part v of the quaternion*/
    const Coordinate& V() const { return _v; } 
    /** @return scalar component w of the quaternion*/
    float& W() { return _w; }
    /** @return scalar component w of the quaternion*/
    const float& W() const { return _w; }
    /** @return the norm of the quaternion*/
    float Norm() const 
    {
      return std::sqrt(_w*_w+_v.x()*_v.x()+_v.y()*_v.y()+_v.z()*_v.z());
    }
    /** invert the quaternion*/
    Quaternion Invert() const
    {
      Quaternion q1;
      
      q1._w=this->_w;
      q1._v.x()=-(this->_v.x());
      q1._v.y()=-(this->_v.y());
      q1._v.z()=-(this->_v.z());
      return q1;
    }
    /** quaternion-quaternion vectorial product*/
    Quaternion operator^(const Quaternion& q)
    {
      Quaternion  q1;
      q1._w=  this->_w*q._w -  this->_v.x()*q._v.x() - this->_v.y()*q._v.y() -this->_v.z()*q._v.z();
      q1._v.x()=this->_w*q._v.x()+this->_v.x()*q._w+this->_v.y()*q._v.z()-this->_v.z()*q._v.y();
      q1._v.y()=this->_w*q._v.y()+this->_v.y()*q._w+this->_v.z()*q._v.x()-this->_v.x()*q._v.z();
      q1._v.z()=this->_w*q._v.z()+this->_v.z()*q._w+this->_v.x()*q._v.y()-this->_v.y()*q._v.x();

      return q1;
    }
    
	/** quaternion-vector product*/
	Quaternion operator^(const Coordinate& c)
	{
	  Quaternion c1(0,c);
	  return (*this)^c1;
	}
	
    void Normalize()
    {
      float norm=this->Norm();
      _v.x()/=norm;
      _v.y()/=norm;
      _v.z()/=norm;
      _w/=norm;
    }
    static Quaternion FromRotationAxis(float angle,const Coordinate& rv)
    {
        angle *= 0.5f;
        Coordinate nrv=rv;
        nrv.Normalize();
 
        float seno=sin(angle);
        Quaternion q;
        q.V().x() = (nrv.x() * seno);
        q.V().y() = (nrv.y() * seno);
        q.V().z() = (nrv.z() * seno);
        q.W() = cos(angle);
       
        return q;
     
    }
    /** return the opengl rotation matrix*/
    void GLMatrix(D2Array<float>& glmatrix)
    {

    float xx      = V().x()*V().x();
    float xy      = V().x() * V().y();
    float xz      = V().x()*V().z();
    float xw      = V().x()*W();
    float yy      = V().y()*V().y();
    float yz      = V().y() * V().z();
    float yw      = V().y()*W();
    float zz      = V().z()*V().z();
    float zw      = V().z() * W();
    glmatrix[0]  = 1 - 2 * ( yy + zz );
    glmatrix[1]  =     2 * ( xy - zw );
    glmatrix[2]  =     2 * ( xz + yw );
    glmatrix[4]  =     2 * ( xy + zw );
    glmatrix[5]  = 1 - 2 * ( xx + zz );
    glmatrix[6]  =     2 * ( yz - xw );
    glmatrix[8]  =     2 * ( xz - yw );
    glmatrix[9]  =     2 * ( yz + xw );
    glmatrix[10] = 1 - 2 * ( xx + yy );
    glmatrix[3]  =glmatrix[7] = glmatrix[11] = glmatrix[12] = glmatrix[13] = glmatrix[14] = 0;
    glmatrix[15] = 1;

    }

private:
    float _w;
    Coordinate _v;

  };
/** used for debugging purposes */

inline KRYOMOLCORE_API std::ostream& operator << (std::ostream& s,const Coordinate& c)
{
  s << "(" << c.x() << "," << c.y() << ","  << c.z() << ")";

  return s;
}

inline KRYOMOLCORE_API std::ostream& operator << ( std::ostream& s, const Quaternion& c)
{
  s << "[" << ( c.W() ) << "," << ( c.V() ) << "]" <<  std::endl;
  return s;
}
}
#endif // KRYOMOLCOORDINATE_H
