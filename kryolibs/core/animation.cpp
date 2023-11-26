/*****************************************************************************************
                            animation.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <qtimer.h>

#include "animation.h"

using namespace kryomol;

Animation::Animation(Molecule& molecule) : m_molecule(molecule)
{

  m_timer = new QTimer(this);
  connect(m_timer,SIGNAL(timeout()),this,SLOT(OnAnimationTimer()));
  m_interval=100;
  m_bfirst=true;
  m_animationmode=0; //first mode in the beginning

}

Animation::~Animation()
{
  m_timer->stop();
}

void Animation::OnAnimationTimer()
{
  if(m_counter>10)
  {
    m_direction*=-1;
    m_counter=0;
  }
  for (size_t i=0;i<m_molecule.CurrentFrame().XYZ().size();++i)
  {
      m_molecule.CurrentFrame().XYZ().at(i) += m_molecule.GetMode().at(m_animationmode).at(i)*m_direction*0.1;
  }

  emit shot();
  m_counter++;
}

void Animation::SetMode(int mode,int frame)
{
  bool brestart=false;
  if( m_timer->isActive() )
  {
    Stop();
    brestart=true;
  }
  m_animationmode=mode;
  if(brestart) Start();
}


void Animation::SetFrame(int frame)
{

  int nframes=frame-m_previousframe[m_animationmode];
  m_previousframe[m_animationmode]=frame;
  const double framescaling=0.10;
  for (size_t i=0;i<m_molecule.CurrentFrame().XYZ().size();++i)
  {
      m_molecule.CurrentFrame().XYZ().at(i) += m_molecule.GetMode().at(m_animationmode).at(i)*nframes*framescaling;
  }
}

void Animation::Start()
{
  if( m_animationmode<0)
    return;
  if(!m_bfirst) Stop();

  m_counter=5;
  m_direction=1;
  m_timer->start(m_interval);
  m_bfirst=false;
}

void Animation::Stop()
{

  m_timer->stop();

  //Fill again the molecule with backup
  RestoreMolecule();
}

void Animation::BackupMolecule()
{
  m_backup.resize(m_molecule.Atoms().size());
  std::vector<Coordinate>::iterator ct=m_backup.begin();
  for(size_t i=0;i<m_molecule.CurrentFrame().XYZ().size();i++,ct++)
  {
    (*ct)=m_molecule.CurrentFrame().XYZ().at(i);
  }

  m_previousframe.resize(m_molecule.Atoms().size());
  std::vector<int>::iterator it;
  for(it=m_previousframe.begin();it!=m_previousframe.end();it++) (*it)=0;

}

void Animation::Break()
{
  m_timer->stop();
  m_bfirst=true;
}

void Animation::RestoreMolecule()
{
  std::vector<Coordinate>::iterator ct=m_backup.begin();
  if (ct==m_backup.end())
    return;
  for(size_t i=0;i<m_molecule.Atoms().size();i++,ct++)
  {
    m_molecule.CurrentFrame().XYZ().at(i)=(*ct);
  }
  std::vector<int>::iterator it;
  for(it=m_previousframe.begin();it!=m_previousframe.end();it++) (*it)=0;
}

bool Animation::IsActive() const
 {
 return m_timer->isActive();
 }

