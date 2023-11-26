#include <sstream>



#include "kryovisoroptical.h"
#include "world.h"

#ifdef Q_OS_MAC
#include <OpenGl/glu.h>
#else
#include <GL/glu.h>
#endif

using namespace kryomol;

KryoVisorOptical::KryoVisorOptical(World* world,QWidget* parent, const QGLWidget* shareWidget, Qt::WindowFlags f ) :
    KryoVisor(world,parent,shareWidget,f), m_bshowelectricdipole(false), m_bshowmagneticdipole(false), m_bshowvelocitydipole(false)
{
    m_activetransition=0;
    m_bshowtransitionchanges=false;
}

void KryoVisorOptical::RenderMolecule(size_t index,GLenum mode)
{
    glPushMatrix();
    Handlers() [index].ApplyTransformation(m_world->Molecules().at(index).CurrentFrameIndex());
    KryoVisor::RenderMolecule(index,mode);
    if ( mode == GL_RENDER)
    {

        Frame& frame=m_world->Molecules().at(index).CurrentFrame();
        if ( m_bshowelectricdipole )
        {
            const kryomol::Coordinate& dipole=frame.GetSpectralLines().at(m_activetransition).ElectricDipole();
            RenderDipole(dipole,QColor(Qt::red));
        }

        if ( m_bshowmagneticdipole )
        {
            const kryomol::Coordinate& dipole=frame.GetSpectralLines().at(m_activetransition).MagneticDipole();
            RenderDipole(dipole,QColor(Qt::green));
        }

        if ( m_bshowvelocitydipole )
        {
            const kryomol::Coordinate& dipole=frame.GetSpectralLines().at(m_activetransition).VelocityDipole();
            RenderDipole(dipole,QColor(255,255,0));
        }
    }
    glPopMatrix();
}

void KryoVisorOptical::OnShowElectricDipole(bool b)
{
    m_bshowelectricdipole=b;
    update();
}

void KryoVisorOptical::OnShowMagneticDipole(bool b)
{
    m_bshowmagneticdipole=b;
    update();
}

void KryoVisorOptical::OnShowVelocityDipole(bool b)
{
    m_bshowvelocitydipole=b;
    update();
}

void KryoVisorOptical::RenderDipole(const kryomol::Coordinate& dipole, const QColor& qcolor)
{

    const int spheres=16;
    GLUquadricObj* quadric= gluNewQuadric();
    static float color[4];
    color[0]=qcolor.redF();
    color[1]=qcolor.greenF();
    color[2]=qcolor.blueF();
    color[3]=qcolor.alphaF();
    glMaterialfv ( GL_FRONT_AND_BACK, GL_AMBIENT, color );
    glMaterialfv ( GL_FRONT_AND_BACK, GL_DIFFUSE, color );
    glPushMatrix();
    {
        Coordinate zaxis ( 0,0,1 );
        Coordinate normal=zaxis^dipole;
        float angle=Coordinate::Angle( zaxis,dipole );
        glRotatef ( angle*180/M_PI,normal.x(),normal.y(),normal.z() );
        float height=dipole.Norm();
        gluCylinder ( quadric,0.05,0.05,height,spheres,spheres );
        glPushMatrix();
        {
            glTranslatef ( 0,0,height );
            gluCylinder ( quadric,0.09,0.001,0.2,spheres,spheres );

        }
        glPopMatrix();

        glPushMatrix();
        {
            glTranslatef ( 0,0,height/2 );
            glDisable ( GL_LIGHTING );

            qglColor ( qcolor );

            //write two spaces to leave some separation between
            //number and atom
            QString sdipole;
            sdipole.sprintf ( "  %.3f a.u.",dipole.Norm());


            RenderText(0.0,0.0,0.0,sdipole);//RenderText ( dipole,0,0, 0 );
            glEnable ( GL_LIGHTING );
        }
        glPopMatrix();
    }
    glPopMatrix();
    gluDeleteQuadric(quadric);

}

void KryoVisorOptical::RenderScreenText()
{
    KryoVisor::RenderScreenText();
    std::stringstream label;
    label << "Current active transition: #" << m_activetransition +1;
    QString str ( label.str().c_str() );
    QColor col ( QColor ( 255,255,255 ) );
    qglColor ( col );
    renderText ( rect().right() -200,rect().bottom() - 30,str,GLFont() );
}

