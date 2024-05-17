/*****************************************************************************************
                            world.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <sstream>

#include "world.h"
#include "plugin.h"
#include "glvisor.h"
#include "kryovisor.h"
#include "kryovisoroptical.h"
//Added by qt3to4:
#include <QDropEvent>
#include <QMouseEvent>
#include <QDragEnterEvent>



using namespace kryomol;

/**  \brief KryoMol simulation world constructor

      A OpenGl GLVisor object will be constructed. This visor will be a child of the QWidget \a parent
      The simulation world is accessed inside the plugins using the GetWorld() method. For instance
      \code
      MyPlugin::Initialize()
      {
        //will print actual frame coodinates in the screen
        std::cout << GetWorld()->CurrentMolecule()->CurrentFrame();
      }
      \endcode
*/
World::World ( QWidget* parent, VisorType vtype, const QGLWidget* shareWidget, Qt::WindowFlags f )
{
    m_currentplugin=NULL;
    m_currentmolecule=0;
    switch (vtype)
    {
        case glvisor:
            m_visor = new GLVisor ( this,parent,shareWidget,f );
        break;
        case optvisor:
            m_visor = new KryoVisorOpt ( this,parent,shareWidget,f  );
        break;
        case freqvisor:
            m_visor = new KryoVisorFreq ( this,parent,shareWidget,f  );
        break;
        case kryovisor:
            m_visor = new KryoVisor ( this,parent,shareWidget,f );
        break;
        case opticalvisor:
            m_visor = new KryoVisorOptical(this,parent,shareWidget,f);
            break;
        default:
            m_visor = new GLVisor ( this,parent,shareWidget,f );
        break;
    }
    m_visor->SelectMolecule ( m_currentmolecule );
}

/**  \brief QryoMol simulation world constructor

      If \a bGUI is set to true a parentless
     OpenGL visor will be also be built
*/
World::World ( bool bGUI ) 
{
  ;
  m_currentplugin=NULL;
  m_currentmolecule=0;

  if ( bGUI )
  {
    m_visor = new GLVisor ( this );
    m_visor->SelectMolecule ( m_currentmolecule );
  }
  else
  {
    m_visor=NULL;
  }

}

World::~World()
{}

/** \brief world initialization

  If a OpenGL visor is associated to the world, this will be also initialized
  */
void World::Initialize()
{
  for(std::vector<Molecule>::iterator mt=m_molecules.begin();mt!=m_molecules.end();++mt)
  {
      SetPopulations(*mt);
      mt->SetColors();
  }

  if ( m_visor ) m_visor->Initialize();
}

/** \brief Set the temperature of the simulation world

   If a NVT ensamble is selected boltzmann populations will be computed accordingly
*/
void World::SetTemperature ( double t )
{
  Thermostat::SetTemperature ( t );
  emit thermostatChanged();
}
/** \brief Set the type of ensamble

  Set the type of ensamble to NVT or NVE type
  In the NVT (constant temperature) case populations will
  be computed accordingly to the boltzmann expression. In the NVE ensamble
  all conformers/frames will have the same population
  */
void World::SetEnsamble ( ensamble e )
{
  Thermostat::SetEnsamble ( e );
  emit thermostatChanged();
}


/** \brief Clear the simulation world

Delete all molecules from the world
*/
void World::Clear()
{
  m_molecules.clear();
}

/** \return A const pointer to the molecule currently active*/
const Molecule* World::CurrentMolecule() const
{
  if ( m_molecules.empty() )
    return NULL;
  else
    return &m_molecules.at ( m_currentmolecule );
}

/** \return A pointer to the molecule currently active*/
Molecule* World::CurrentMolecule()
{
  if ( m_molecules.empty() )
    return NULL;
  else
    return &m_molecules.at ( m_currentmolecule );
}
/** \return the index (0 to N-1 ) of the molecule currently active*/
size_t World::CurrentMoleculeIndex() const
{
  return m_currentmolecule;
}

/** \brief select a frame for the current molecule

   Select conformer \a frame for the molecule currently active
   This method will emit the currentFrame(size_t ) signal
   */
void  World::SelectFrame ( size_t frame )
{
  if (!CurrentMolecule() ) return;

    if ( frame >= CurrentMolecule()->Frames().size() )
    {
      std::cerr << "World:: Invalid frame index" << frame << std::endl;
      return;
    }

    CurrentMolecule()->SetCurrentFrame ( frame );

    if ( m_visor ) m_visor->Center();

    emit currentFrame ( frame );

}

/** \brief select a molecule inside the world*

  Select the molecule of index \a mol
  This method will emit the signal currentMolecule(size_t )
  */
void  World::SelectMolecule ( size_t mol )
{
  if ( mol >= m_molecules.size() )
  {
    std::cerr << "World:: Invalid molecule index" << mol << std::endl;
    return;
  }
  m_currentmolecule=mol;


  if ( m_visor ) m_visor->update();

  emit currentMolecule ( mol );
}

/** \return a reference to a stl vector of Molecule objects */
std::vector<Molecule>& World::Molecules()
{
  return m_molecules;
}

/** \return a reference to a stl vector of all Molecule objects */
const std::vector<Molecule>& World::Molecules() const
{
  return m_molecules;
}

/** \return a pointer to the associated GLVisor */
GLVisor* World::Visor()
{
  return m_visor;
}

/** \return a const pointer to the associated GLVisor */
const GLVisor* World::Visor() const
{
  return m_visor;
}


/** Associate the current active plugin \a p with the simulation world*/
void World::SetCurrentPlugin ( Plugin* p )
{
  m_currentplugin=p;
}

/** \return a pointer to the current plugin associated to the simulation world*/
Plugin* World::CurrentPlugin()
{
  return m_currentplugin;
}

/** \return a const pointer to the plugin associated to the simulation world*/
const Plugin* World::CurrentPlugin() const
{
  return m_currentplugin;
}
