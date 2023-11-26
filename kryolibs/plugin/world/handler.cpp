/*****************************************************************************************
                            handler.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <QtOpenGL/QGLWidget>
#include "handler.h"


namespace kryomol
{
class FrameHandler
{
  public:
    FrameHandler () :  m_binitial ( true )
    {
      m_matrix=D2Array<float> ( 4,4,0 );
      m_matrix ( 0,0 ) =m_matrix ( 1,1 ) =m_matrix ( 2,2 ) =m_matrix ( 3,3 ) =1;
      m_quat=Quaternion ( 0,Coordinate ( 1,1,1 ) );
    }
    ~FrameHandler() {}
    float m_scale;
    Coordinate m_rotcenter;
    D2Array<float> m_matrix;
    Coordinate m_rotvector;
    Quaternion m_quat;
    Quaternion m_rotangle;
    bool m_binitial;
};
}

using namespace kryomol;

class kryomol::MoleculeHandlerPrivate
{
  public:
  MoleculeHandlerPrivate() : m_scale(0.3), m_blocked(true) {}
  ~MoleculeHandlerPrivate() {}
  std::vector<FrameHandler> m_handlers;
  float m_scale;
  bool m_blocked;
};

/** \return a const reference to the 4x4 rotation matrix for \a frame*/
const D2Array<float>& MoleculeHandler::Matrix(size_t frame) const
{
  return _d->m_handlers.at(frame).m_matrix;
}

/** \return a reference to the 4x4 rotation matrix for \a frame*/
D2Array<float>& MoleculeHandler::Matrix(size_t frame)
{
  return _d->m_handlers.at(frame).m_matrix;
}

/** rotate the conformer \a frame by the angle \a rotangle (radians) around the vector \a rotvector*/
void MoleculeHandler::RotateFrame ( float rotangle,const Coordinate& rotvector,size_t frame )
{
  rotangle*=M_PI/180.;
  Quaternion q=Quaternion::FromRotationAxis ( rotangle,rotvector );
  q.Normalize();

  if ( _d->m_handlers[frame].m_binitial )
  {
    _d->m_handlers[frame].m_quat=q;
    _d->m_handlers[frame].m_binitial=false;
  }
  else
  {
    _d->m_handlers[frame].m_quat=_d->m_handlers[frame].m_quat^q;
  }
  _d->m_handlers[frame].m_quat.Normalize();
  _d->m_handlers[frame].m_quat.GLMatrix ( _d->m_handlers[frame].m_matrix );

}

/** rotate the conformer \a frame by the angle \a rotangle (radians) around the vector \a rotvector
if handler is blocked all the conformers will be rotated*/
void MoleculeHandler::Rotate ( float rotangle,const Coordinate& rotvector,size_t frame )
{
  if ( _d->m_blocked )
  {
      for(size_t i=0;i<_d->m_handlers.size();++i)
          RotateFrame(rotangle,rotvector,i);
  }
  else
  {
      RotateFrame( rotangle,rotvector,frame);
  }
}

/** Apply transformation matrix to the conformer \a frame*/
void MoleculeHandler::ApplyTransformation(size_t frame)
{

  glLoadIdentity();
  FrameHandler& f=_d->m_handlers[frame];
  glTranslatef ( f.m_rotcenter.x(),f.m_rotcenter.y(),f.m_rotcenter.z() );

  glMultMatrixf ( ( GLfloat* ) f.m_matrix );
  glScalef ( _d->m_scale, _d->m_scale, _d->m_scale );
  glTranslatef ( -f.m_rotcenter.x(),-f.m_rotcenter.y(),-f.m_rotcenter.z() );

}

/** Set the trackball rotation center for conformer \a frame or globally if blocked*/
void MoleculeHandler::SetRotationCenter( const Coordinate& c,size_t frame )
{
   /* if ( _d->m_blocked )
    {
        for(std::vector<FrameHandler>::iterator it=_d->m_handlers.begin();it!=_d->m_handlers.end();++it)
            it->m_rotcenter=c;
    }
    else*/
     _d->m_handlers[frame].m_rotcenter=c;
}


/** \return a reference ro the trackball rotation center*/
const Coordinate& MoleculeHandler::RotationCenter(size_t frame) const
{
  return _d->m_handlers[frame].m_rotcenter;
}

/** built a MoleculeHandler object*/
MoleculeHandler::MoleculeHandler()
{
  _d = new MoleculeHandlerPrivate();
}

/** \brief copy constructor*/
MoleculeHandler::MoleculeHandler(const MoleculeHandler& other)
{
  _d= new MoleculeHandlerPrivate();
}

MoleculeHandler::~MoleculeHandler()
{
  delete _d;
  
}

/** \brief copy operator*/
MoleculeHandler& MoleculeHandler::operator=(const MoleculeHandler& other)
{
  //do absolutely nothing
  return *this; 
}

/** \return A reference to the scale*/
float& MoleculeHandler::Scale()
{
  return _d->m_scale;
} 

/** \return A const referene to the scale*/
const float& MoleculeHandler::Scale() const
{
  return _d->m_scale;
}

/** \brief Set the number of frames/conformers*/
void MoleculeHandler::SetNFrames(size_t nframes)
{
  _d->m_handlers.resize(nframes);
}

/** \return the number of frames/conformers associated to the handler*/
size_t MoleculeHandler::NFrames() const
{
  return _d->m_handlers.size();
}

/** \return true if rotation handler is blocked to the same value for all conformersr*/
bool MoleculeHandler::IsBlocked() const
{
  return _d->m_blocked;
}

/** set a common rotation state for all conformers*/
void MoleculeHandler::SetBlocked(bool b)
{
  _d->m_blocked=b;
}

/** reset the trackball rotation state*/
void MoleculeHandler::Reset()
{
    for(std::vector<FrameHandler>::iterator it=_d->m_handlers.begin();it!=_d->m_handlers.end();++it)
    {
        (*it)=FrameHandler();
    }
}
