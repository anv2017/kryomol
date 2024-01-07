/*****************************************************************************************
                            glvisor.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <QToolBar>
#include <sstream>
#include <algorithm>
#include <QActionGroup>

#ifdef Q_OS_MAC
#include <OpenGl/glu.h>
#else
#include <GL/glu.h>
#endif

#include "glvisor.h"
#include "gl2ps.h"
#include "stringtools.h"
#include "mathtools.h"
#include "orbitalarray.h"
#include "world.h"
#include "plugin.h"
#include "qryomolapp.h"
#include "density.h"
#include "renderdensity.h"

using namespace kryomol;

const float rotationpass=5.0f;

class GLVisor::GLVisorPrivate
{
public:
    GLVisorPrivate() {}
    ~GLVisorPrivate() {}
    //draw wireframe when moving
    bool m_bwfonmoving;
};

/** \brief Constructor
 Construct a visor associated to the World world \a world as a child of the QWidget \a parent
*/
GLVisor::GLVisor ( World* world,QWidget* parent, const QGLWidget* shareWidget, Qt::WindowFlags f )
    :GLVisorBase ( parent,shareWidget,f ) , m_world ( world ) , m_selmode ( NONE )
{

    _d = new GLVisorPrivate();
    _d->m_bwfonmoving=true;
    m_bdistancechange = false;
    m_banglechange = false;
    m_bdihedralchange = false;
    m_bshowdistances = false;
    m_bshowdipole = false;
    m_bshowcell = false;
    m_bshowaxis = false;
    m_bshowdensity = false;
    m_transparence = 0.5;

    setAcceptDrops ( true );
    InitToolBars();
}

/** \brief Destructor
*/
GLVisor::~GLVisor()
{
    delete _d;
}

/** \brief Init GLVisor associated toolbars

A toolbar with actions such as molecular measurements.. is constructed here
*/
void GLVisor::InitToolBars()
{

    QToolBar* measurebar= new QToolBar ();
    QActionGroup* measureActions = new QActionGroup ( this );
    QAction* noneAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/select.png") ),tr ( "Plain Selection" ),this );
    noneAction->setActionGroup ( measureActions );
    connect ( noneAction,SIGNAL ( activated() ),this,SLOT ( OnNoneSelection() ) );
    QAction* distanceAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/distance.png") ),tr ( "Measure Distances" ),this );
    distanceAction->setActionGroup ( measureActions );
    connect ( distanceAction,SIGNAL ( activated() ),this,SLOT ( OnMeasureDistances () ) );

    QAction* angleAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/angle.png") ), tr ( "Measure Angles" ),this );
    connect ( angleAction,SIGNAL ( activated() ),this,SLOT ( OnMeasureAngles () ) );
    angleAction->setActionGroup ( measureActions );
    QAction* dihedralAction = new QAction ( QIcon ( QString::fromUtf8(":/icons/dihedral.png") ),tr ( "Measure Dihedrals" ),this );
    dihedralAction->setActionGroup ( measureActions );
    connect ( dihedralAction,SIGNAL ( activated() ),this,SLOT ( OnMeasureDihedrals () ) );

    noneAction->setCheckable ( true );
    distanceAction->setCheckable ( true );
    angleAction->setCheckable ( true );
    dihedralAction->setCheckable ( true );
    measureActions->setExclusive ( true );


    measurebar->addAction ( noneAction );
    measurebar->addAction ( distanceAction );
    measurebar->addAction ( angleAction );
    measurebar->addAction ( dihedralAction );
    noneAction->setChecked ( true );

    m_toolbars.push_back ( measurebar );

}

/** \brief GLVisor initialization

Initialize the visor accordingly with molecular information
*/
void GLVisor::Initialize()
{
    if ( m_world->CurrentMolecule()->Atoms().size() > 100 )
    {

        _d->m_bwfonmoving=true;
    }
    else
        _d->m_bwfonmoving=false;

    Handlers().resize( m_world->Molecules().size() );

    for(size_t i=0;i<m_world->Molecules().size();++i)
    {
        Handlers()[i].SetNFrames(m_world->Molecules()[i].Frames().size());
    }
    CenterMolecule ( true,true );
    SetGraphMode ( WIREFRAME );
    update();
}

/** manage drag events*/
void GLVisor::dragEnterEvent ( QDragEnterEvent* e )
{
    e->accept();
}

/** manage drops events*/
void GLVisor::dropEvent ( QDropEvent* e )
{
    emit sourceDropped ( e->mimeData() );
}

/** \brief Render OpenGL scene

This method will call first the GLVisorBase::RenderScene methid
and then will call the RenderPlugins() and RenderScreenText() methods
*/
void GLVisor::RenderScene ( GLenum mode )
{
    GLVisorBase::RenderScene ( mode );

    glPushMatrix();

    for ( size_t i=0;i<m_world->Molecules().size();++i )
    {
        glPushMatrix();
        RenderMolecule ( i,mode );
        glPopMatrix();
    }

    glPopMatrix();

    if ( m_world->CurrentMolecule() )
    {
        glPushMatrix();
        Handlers() [m_world->CurrentMoleculeIndex()].ApplyTransformation(m_world->CurrentMolecule()->CurrentFrameIndex());
        RenderPlugins ( mode );
        glPopMatrix();
        if ( mode == GL_RENDER )
        {
            RenderScreenText();
        }
    }

}

/** \brief Do custom rendering for active plugin

*/
void GLVisor::RenderPlugins ( GLenum mode )
{

    static float icolor[4];
    icolor[0]=0.9f;
    icolor[1]=0.3f;
    icolor[2]=0.3f;
    icolor[3]=0.5f;

    if ( mode == GL_RENDER && m_selmode != NONE )
    {
        const Molecule* mol=m_world->CurrentMolecule();
        if ( mol )
        {
            GLUquadricObj* quadric= gluNewQuadric();
            glMaterialfv ( GL_FRONT_AND_BACK, GL_AMBIENT, icolor );
            glMaterialfv ( GL_FRONT_AND_BACK, GL_DIFFUSE, icolor );

            for ( size_t i=0;i<m_selatoms.size();++i )
            {
                glPushMatrix();
                const Coordinate& c=mol->CurrentFrame().XYZ() [m_selatoms[i]];
                glTranslatef ( c.x(),c.y(),c.z() );
                gluSphere ( quadric,0.25,8,8 );
                glPopMatrix();
            }
            gluDeleteQuadric ( quadric );
        }
    }
    //Ok load the plugin
    if ( mode == GL_RENDER )
    {

        if ( m_world->CurrentPlugin() )
        {
            m_world->CurrentPlugin()->Render();
        }
    }
}

/** \brief Render text on the opengl screen

 This function will render the current index of the selected conformer
*/
void GLVisor::RenderScreenText()
{
    if ( !m_world->CurrentMolecule() ) return;
    std::stringstream label;
    label << "#" << m_world->CurrentMolecule()->CurrentFrameIndex() +1 << " of " << m_world->CurrentMolecule()->Frames().size();
    QString str ( label.str().c_str() );

    int left=this->rect().left();
    int bottom=this->rect().bottom();

    qglColor ( QColor ( 255,255,255 ) );
    RenderText ( left+40,bottom-40, str );
}

/** Process OpenGL selection events

This function will detect hits over rendered atoms and will call
the Plugin::ProcessSelection(int ) method for the current plugin
*/
void GLVisor::OnSelection ( QMouseEvent* e )
{

    static GLuint selectBuff[64];
    GLint hits, viewport[4];
    glSelectBuffer ( 64,selectBuff );
    glGetIntegerv ( GL_VIEWPORT,viewport );
    glMatrixMode ( GL_PROJECTION );
    glPushMatrix();
    glRenderMode ( GL_SELECT );
    glLoadIdentity();

    //there is a problem at least with some macosx systems and Qt 5.3
    //the viewport is scaled as twice the height() of the system.
    //so i have to detect at run time and scale everything accordingly
    int scale=viewport[3]/this->height();
    //qDebug() << "scale=" << scale;
    //An 8x8 pick seems to be necessary on win and macosx
    gluPickMatrix ( scale*e->x(),viewport[3]-scale*e->y(),scale*8,scale*8,viewport );
    SetupPerspective();
    glMatrixMode ( GL_MODELVIEW );
    //We have to multiply here by the rotation acumulated matrix
    //glMultMatrixf ( ( GLfloat * ) m_transform );

    RenderScene ( GL_SELECT );
    hits=glRenderMode ( GL_RENDER );
    //OK lets call the current plugin and process the selection
    if ( hits > 0 ) //dont use (hits) since it can be -1
    {
        int i;
        GLuint names, *ptr, minZ,*ptrNames, numberOfNames;

        ptr = ( GLuint * ) selectBuff;
        minZ = 0xffffffff;
        for ( i = 0; i < hits; i++ )
        {
            names = *ptr;
            ptr++;
            if ( *ptr < minZ )
            {
                numberOfNames = names;
                minZ = *ptr;
                ptrNames = ptr+2;
            }

            ptr += names+2;
        }

        ptr = ptrNames;

        int natom=*ptr;
        if ( m_selmode == NONE )
        {
            if ( m_world->CurrentPlugin() )
                m_world->CurrentPlugin()->ProcessSelection ( natom );
        }
        else
            ProcessOwnSelection ( natom );
    }

    glMatrixMode ( GL_PROJECTION );
    glPopMatrix();

    glMatrixMode ( GL_MODELVIEW );

    if ( hits > 0 )
        update();

}

/** \brief vector picture exporting

Export a vectore picture to the file \a file using the format \a mode
The exporting is made through the gl2ps library. Current available formats are
 PS, EPS,SVG, PDF
*/
void GLVisor::ExportVectorPicture ( const QString& file, const QString& mode )
{
    std::cout << "Exporting in " << mode.toStdString() << " format to file " << file.toStdString() << std::endl;
    QString upmode=mode.toUpper();
    char *oldlocale = setlocale ( LC_NUMERIC, "C" );
    FILE* fp=fopen ( file.toStdString().c_str(),"wb" );
    if ( !fp ) return;
    GLint buffsize=0, state =GL2PS_OVERFLOW;
    GLint viewport[4];
    viewport[0]=0;
    viewport[1]=0;
    viewport[2]=640;
    viewport[3]=480;
    glGetIntegerv ( GL_VIEWPORT, viewport );

    int gmode=GL2PS_PS;

    if ( upmode == "EPS" ) gmode=GL2PS_EPS;
    if ( upmode =="PS" ) gmode =GL2PS_PS;
    if ( upmode == "SVG" ) gmode=GL2PS_SVG;
    if ( upmode == "PDF" ) gmode =GL2PS_PDF;

    SetVectorGraphicsMode ( true );
    while ( state == GL2PS_OVERFLOW )
    {

        buffsize+=1024*1024;
        gl2psBeginPage ( file.toStdString().c_str(), "QryoMol", viewport,
                         gmode , GL2PS_BSP_SORT, GL2PS_DRAW_BACKGROUND | GL2PS_SILENT | GL2PS_BEST_ROOT ,
                         GL_RGBA, 0, NULL, 0, 0, 0, buffsize,
                         fp,file.toStdString().c_str() );
        gl2psEnable ( GL2PS_BLEND );
        RenderScene();
        gl2psDisable ( GL2PS_BLEND );
        state=gl2psEndPage();

    }
    SetVectorGraphicsMode ( false );
    fclose ( fp );
    setlocale ( LC_NUMERIC, oldlocale );

}

QSize GLVisor::sizeHint()
{
    return QSize ( 800,600 );
}


/** handle mouse press events*/
void GLVisor::mousePressEvent ( QMouseEvent* e )
{
    GLVisorBase::mousePressEvent ( e );
    switch ( e->button() )
    {
    case Qt::LeftButton:
        OnSelection ( e );
        break;
    case Qt::RightButton:
        break;
    default:
        break;
    }
}

/** handle mouse move events*/
void GLVisor::mouseMoveEvent ( QMouseEvent* e )
{
    if ( !m_world->CurrentMolecule() ) return;

    int x=MousePosition().x()-e->pos().x();
    int y=MousePosition().y()-e->pos().y();

    if ( abs ( x ) <2 && abs ( y ) <2 )
        return;

    if ( GrabbedButtons().testFlag ( Qt::LeftButton ) &&
         !GrabbedButtons().testFlag ( Qt::RightButton ) )
    {
        Coordinate c ( x,-y,0 );
        //Get a perpendicular vector
        Coordinate rotationvector=c^Coordinate ( 0,0,1 );

#ifdef __GNUC__
#warning multimolecule
#endif
        if ( !Handlers().empty() )
        {
            Handlers() [m_world->CurrentMoleculeIndex() ].Rotate ( -rotationpass,rotationvector,m_world->CurrentMolecule()->CurrentFrameIndex() );
            updateGL();
        }

    }
    if ( GrabbedButtons().testFlag ( Qt::RightButton ) &&
         !GrabbedButtons().testFlag ( Qt::LeftButton ) )
    {
        if ( MousePosition().y() < e->pos().y() )
        {

            Handlers() [m_world->CurrentMoleculeIndex()].Scale() *=1.3;

        }
        if ( MousePosition().y() > e->pos().y() )
        {
            Handlers() [m_world->CurrentMoleculeIndex()].Scale() /=1.3;
        }
        update();
    }

    if ( GrabbedButtons().testFlag ( Qt::MidButton )
         || ( GrabbedButtons().testFlag ( Qt::LeftButton ) && GrabbedButtons().testFlag ( Qt::RightButton ) ) )
    {
        if ( x != 0 || y != 0 )
        {
            Coordinate newc=Handlers() [m_world->CurrentMoleculeIndex()].RotationCenter(m_world->CurrentMolecule()->CurrentFrameIndex());

            if ( x != 0 ) newc.x() -= ( x*0.1 );
            if ( y != 0 ) newc.y() += ( y*0.1 );


            Handlers() [m_world->CurrentMoleculeIndex()].SetRotationCenter ( newc,m_world->CurrentMolecule()->CurrentFrameIndex() );
            update();
        }
    }

    GLVisorBase::mouseMoveEvent ( e );

}


/** Set the selection mode to DISTANCES*/
void GLVisor::OnMeasureDistances ()
{
    m_selmode=DISTANCES;
    update();
}

/** Set the selection mode to ANGLES*/
void GLVisor::OnMeasureAngles ()
{
    m_selmode=ANGLES;
    update();
}

/** Set the selection mode to DIHEDRALS*/
void GLVisor::OnMeasureDihedrals ()
{
    m_selmode=DIHEDRALS;
    update();
}

/**Set the selection mode to generic four atoms*/
void GLVisor::OnFourAtomsSelectionMode()
{
    m_selmode=GENERICFOURATOMS;
    update();
}

/** Set the selection mode to ROTATEBOND*/
void GLVisor::OnRotateBond()
{
    m_selmode=ROTATEBOND;
    update();
}


/** Set the selection mode to NONE*/
void GLVisor::OnNoneSelection()
{
    m_selmode=NONE;
    update();
}

/** process selection event if selection mode is not NONE*/
void GLVisor::ProcessOwnSelection ( int atom )
{
    std::vector<size_t>::iterator it=std::find ( m_selatoms.begin(),m_selatoms.end(),static_cast<size_t> ( atom ) );
    if ( it != m_selatoms.end() )
    {
        m_selatoms.erase ( it );
        return;
    }
    m_selatoms.push_back ( atom );
    const Molecule& mol=*m_world->CurrentMolecule();
    const Frame& fr=m_world->CurrentMolecule()->CurrentFrame();
    std::cout << "m_selmode=" << m_selmode << std::endl;
    switch ( m_selmode )
    {
    case DISTANCES:
        if ( m_selatoms.size() == 2 )
        {
            QString str;
            str.sprintf ( "(%s%d,%s%d)  %.3f",
                          mol.Atoms() [m_selatoms[0]].Symbol().c_str(),
                    m_selatoms[0]+1,
                    mol.Atoms() [m_selatoms[1]].Symbol().c_str(),
                    m_selatoms[1]+1,
                    fr.Distance ( m_selatoms[0],m_selatoms[1] ) );
            m_distances.push_back ( Molecule::pair ( m_selatoms[0],m_selatoms[1],str.toStdString() ) );
            emit distance ( str );
            OnResetSelection ();
        }
        break;
    case ANGLES:
        if ( m_selatoms.size() == 3 )
        {
            QString str;
            str.sprintf ( "(%s%d,%s%d,%s%d)  %.2f",
                          mol.Atoms() [m_selatoms[0]].Symbol().c_str(),
                    m_selatoms[0]+1,
                    mol.Atoms() [m_selatoms[1]].Symbol().c_str(),
                    m_selatoms[1]+1,
                    mol.Atoms() [m_selatoms[2]].Symbol().c_str(),
                    m_selatoms[2]+1,
                    180*fr.Angle ( m_selatoms[0],m_selatoms[1],m_selatoms[2] ) /M_PI );
            m_angles.push_back ( Molecule::triad ( m_selatoms[0],m_selatoms[1],m_selatoms[2],str.toStdString() ) );
            //emit text ( str );
            emit angle ( str );
            OnResetSelection ();
        }
        break;
    case DIHEDRALS:
        if ( m_selatoms.size() == 4 )
        {
            QString str;
            str.sprintf ( "(%s%d,%s%d,%s%d,%s%d)  %.2f",
                          mol.Atoms() [m_selatoms[0]].Symbol().c_str(),
                    m_selatoms[0]+1,
                    mol.Atoms() [m_selatoms[1]].Symbol().c_str(),
                    m_selatoms[1]+1,
                    mol.Atoms() [m_selatoms[2]].Symbol().c_str(),
                    m_selatoms[2]+1,
                    mol.Atoms() [m_selatoms[3]].Symbol().c_str(),
                    m_selatoms[3]+1,
                    180*fr.Dihedral ( m_selatoms[0],m_selatoms[1],m_selatoms[2],m_selatoms[3] ) /M_PI );
            m_dihedrals.push_back ( Molecule::quad ( m_selatoms[0],m_selatoms[1],m_selatoms[2],
                    m_selatoms[3],str.toStdString() ) );

            //emit text ( str );
            emit dihedral ( str );
            OnResetSelection ();
        }
        break;
    case GENERICFOURATOMS:
        if ( m_selatoms.size() == 4)
        {
            emit selectedatoms(m_selatoms);
            OnResetSelection();
        }
        break;
    case ROTATEBOND:
        if ( m_selatoms.size() == 2 )
        {
            QString str;
            //m_rotbond=( Molecule::pair ( m_selatoms[0],m_selatoms[1],str.toStdString() ) );
            //OnResetSelection ();
        }
        break;
    default:
        break;

    }
}

/** reset the selected atoms*/
void GLVisor::OnResetSelection()
{
    m_selatoms.clear();
    update();
}

/** \brief do molecular rendering

This method does actual rendering of molecules.
Current modes are WIREFRAME, CPK and Sticks
*/


void GLVisor::RenderMolecule ( size_t index,GLenum mode/*=GL_RENDER*/ )
{
    const Molecule& molecule= m_world->Molecules() [index];
    if ( molecule.Frames().empty() ) return;

    glPushMatrix();
    Handlers()[m_world->CurrentMoleculeIndex()].ApplyTransformation(m_world->CurrentMolecule()->CurrentFrameIndex());

    GLUquadricObj* quadric=gluNewQuadric();
    gluQuadricDrawStyle ( quadric, ( GLenum ) GLU_FILL );
    graphmode gmode;

    if ( IsMouseMoving() && _d->m_bwfonmoving )
        gmode=WIREFRAME;
    else gmode=GraphMode();


    glInitNames();
    glPushName ( 0 );
    std::vector<Atom>::const_iterator mit;
    const Frame& frame=molecule.CurrentFrame();

    std::vector<Coordinate>::const_iterator ct=frame.XYZ().begin();
    int i=0;

    if ( ! ( gmode == WIREFRAME && mode == GL_RENDER && !ShowSymbols() && !ShowNumbers() && !ShowPDBInfo()))// && !ShowDipole() && !ShowCell())) //dont draw anything for wireframe in render mode
        for ( mit=molecule.Atoms().begin();mit!=molecule.Atoms().end();++mit,i++,++ct )
        {
            bool  visible=true;
            if ( !molecule.Residues().empty() )
            {
                visible=mit->Residue()->Visible();
            }
            if ( visible )
            {

                glPushMatrix();
                glTranslatef ( ct->x(), ct->y(), ct->z() );
                if ( mode == GL_SELECT )
                    glLoadName ( i );
                GetMaterialColor ( *mit );
                if ( mode == GL_RENDER )
                {
                    if ( ShowNumbers() || ShowSymbols() || ShowPDBInfo())// || ShowDipole() || ShowCell())
                    {
                        if ( !IsMouseMoving() )
                        {
                            glDisable ( GL_LIGHTING );

                            static QColor col ( QColor ( 255,255,255 ) );
                            qglColor ( col );

                            QString number="   ";
                            if ( ShowSymbols() )
                                number+=mit->Symbol().c_str();
                            if ( ShowNumbers() )
                            {
                                QString snum;
                                snum.sprintf ( "%d",i+1 );
                                number+=snum;
                            }

                            if ( ShowPDBInfo() )
                            {
                                if ( mit->Residue() )
                                {
                                    QString pdb=QString(mit->PDBName().c_str()) +","+QString(mit->Residue()->Name().c_str())+QString(mit->Residue()->Index().c_str());
                                    number+=pdb;
                                }
                            }

                            RenderText ( 0,0, 0, number,GLFont() );
                            glEnable ( GL_LIGHTING );
                        }
                    }
                }

                if ( gmode == CPK )
                    gluSphere ( quadric,mit->VdW() /*0.05*/,2*Spheres(),2*Spheres() );
                if ( gmode == STICKS )
                    gluSphere ( quadric,0.1,Spheres(),Spheres() );

                if ( gmode == WIREFRAME && mode == GL_SELECT ) //this is necessary for selection only but slows down a lot the drawing )
                {
                    glBegin ( GL_POINTS );
                    glVertex3f ( 0,0,0 );
                    glEnd();
                }
                glPopMatrix();

            }
        }


    if ( mode == GL_RENDER )
    {
#ifdef __GNUC__
#warning multimolecule
#endif
        glLineWidth ( 1.0*Handlers() [m_world->CurrentMoleculeIndex() ].Scale() );

        if ( gmode == WIREFRAME )
        {
            glDisable ( GL_LIGHTING );
            glBegin ( GL_LINES );
            std::vector<Bond>::const_iterator cit;
            const std::vector<Coordinate>& c=frame.XYZ();

            if (molecule.Atoms().size() == 1)
            {
                GetColor( molecule.Atoms().back() );
                glVertex3f( c.back().x()-0.15, c.back().y(), c.back().z());
                glVertex3f( c.back().x()+0.15, c.back().y(), c.back().z());
                glVertex3f( c.back().x(), c.back().y()-0.15, c.back().z());
                glVertex3f( c.back().x(), c.back().y()+0.15, c.back().z());
            }
            else
            {
                for ( cit=molecule.Bonds().begin();cit!=molecule.Bonds().end();++cit )
                {
                    bool visible=true;
                    if ( !molecule.Residues().empty() )
                        visible= ( molecule.Atoms() [cit->I() ].Residue()->Visible() && molecule.Atoms() [cit->J() ].Residue()->Visible() );

                    if ( visible )
                    {
                        if ( molecule.Atoms() [ cit->I() ].Z() == molecule.Atoms() [ cit->J() ].Z() )
                        {

                            GetColor ( molecule.Atoms() [cit->I() ] );
                            glVertex3f ( c[cit->I() ].x(),c[cit->I() ].y(),c[cit->I() ].z() );
                            glVertex3f ( c[cit->J() ].x(),c[cit->J() ].y(),c[cit->J() ].z() );
                        }
                        else
                        {
                            Coordinate middle=Coordinate::MiddlePoint ( c[cit->I() ],c[cit->J() ] );
                            GetColor ( molecule.Atoms() [cit->I() ] );
                            glVertex3f ( c[cit->I() ].x(),c[cit->I() ].y(),c[cit->I() ].z() );
                            glVertex3f ( middle.x(),middle.y(),middle.z() );
                            GetColor ( molecule.Atoms() [cit->J() ] );
                            glVertex3f ( middle.x(),middle.y(),middle.z() );
                            glVertex3f ( c[cit->J() ].x(),c[cit->J() ].y(),c[cit->J() ].z() );
                        }
                    }
                }
            }

            glEnd();
            glEnable ( GL_LIGHTING );
        }

        if ( gmode == STICKS )
        {
            std::vector<Bond>::const_iterator cit;
            const std::vector<Coordinate>& c=frame.XYZ();
            for ( cit=molecule.Bonds().begin();cit!=molecule.Bonds().end();++cit )
            {
                bool visible=true;
                if ( !molecule.Residues().empty() )
                    visible= ( molecule.Atoms() [cit->I() ].Residue()->Visible() && molecule.Atoms() [cit->J() ].Residue()->Visible() );

                if ( visible )
                {
                    if ( molecule.Atoms() [ cit->I() ].Z() == molecule.Atoms() [ cit->J() ].Z() )
                    {
                        glPushMatrix();
                        glTranslatef ( c[cit->I() ].x(),c[cit->I() ].y(),c[cit->I() ].z() );
                        Coordinate bond= c[cit->J() ]-c[cit->I() ];
                        GetMaterialColor ( molecule.Atoms() [cit->I() ] );
                        Coordinate zaxis ( 0,0,1 );
                        Coordinate normal=zaxis^bond;
                        float angle=Coordinate::Angle ( zaxis,bond );
                        glRotatef ( angle*180/M_PI,normal.x(),normal.y(),normal.z() );
                        gluCylinder ( quadric,0.1,0.1, bond.Norm() , Spheres(),Spheres() );
                        glPopMatrix();
                    }
                    else
                    {
                        glPushMatrix();
                        glTranslatef ( c[cit->I() ].x(),c[cit->I() ].y(),c[cit->I() ].z() );
                        Coordinate bond= c[cit->J() ]-c[cit->I() ];
                        float length=bond.Norm() /2;
                        GetMaterialColor ( molecule.Atoms() [cit->I() ] );
                        Coordinate zaxis ( 0,0,1 );
                        Coordinate normal=zaxis^bond;
                        float angle=Coordinate::Angle ( zaxis,bond );
                        glRotatef ( angle*180/M_PI,normal.x(),normal.y(),normal.z() );
                        gluCylinder ( quadric,0.1,0.1, length , Spheres() , Spheres() );
                        glPopMatrix();

                        glPushMatrix();
                        glTranslatef ( c[cit->J() ].x(),c[cit->J() ].y(),c[cit->J() ].z() );
                        bond= c[cit->I() ]-c[cit->J() ];
                        GetMaterialColor ( molecule.Atoms() [cit->J() ] );
                        normal=zaxis^bond;
                        angle=Coordinate::Angle ( zaxis,bond );
                        glRotatef ( angle*180/M_PI,normal.x(),normal.y(),normal.z() );
                        gluCylinder ( quadric,0.1,0.1, length , Spheres() , Spheres() );
                        glPopMatrix();

                    }
                }
            }
        }
    }

    if (m_bdistancechange)
    {
        Molecule::pair atompair = m_distances [m_pair];

        static float color[4];
        color[0]=0.9f;
        color[1]=0.3f;
        color[2]=0.3f;
        color[3]=0.5f;

        glMaterialfv ( GL_FRONT_AND_BACK, GL_AMBIENT, color );
        glMaterialfv ( GL_FRONT_AND_BACK, GL_DIFFUSE, color );

        glPushMatrix();
        glTranslatef ( molecule.CurrentFrame().XYZ()[atompair.i].x(),
                molecule.CurrentFrame().XYZ()[atompair.i].y(),
                molecule.CurrentFrame().XYZ()[atompair.i].z()  );
        gluSphere(quadric, 0.25, 8, 8);
        glPopMatrix();
        glPushMatrix();
        glTranslatef( molecule.CurrentFrame().XYZ()[atompair.j].x(),
                molecule.CurrentFrame().XYZ()[atompair.j].y(),
                molecule.CurrentFrame().XYZ()[atompair.j].z() );
        gluSphere(quadric, 0.25, 8, 8);
        glPopMatrix();

        if ( m_bshowdistances )
        {

            glLineWidth ( 1 );
            glEnable ( GL_LINE_STIPPLE );
            glLineStipple ( 1,0xCCCC );
            glDisable ( GL_LIGHTING );

            //draw negative values in blue
            if ( atompair.value.at( 0 ) == '-' )
                glColor3f ( 0.5,0.5,0.9 );
            else
                glColor3f ( 0.9,0.5,0.5 );

            glBegin ( GL_LINES );
            glVertex3f ( molecule.CurrentFrame().XYZ()[atompair.i].x(),
                    molecule.CurrentFrame().XYZ()[atompair.i].y(),
                    molecule.CurrentFrame().XYZ()[atompair.i].z()  );

            glVertex3f ( molecule.CurrentFrame().XYZ()[atompair.j].x(),
                    molecule.CurrentFrame().XYZ()[atompair.j].y(),
                    molecule.CurrentFrame().XYZ()[atompair.j].z() );

            Coordinate middle=Coordinate::MiddlePoint ( molecule.CurrentFrame().XYZ()[atompair.i],molecule.CurrentFrame().XYZ()[atompair.j]);
            glEnd();

            //write two spaces to leave some separation between
            //number and atom
            QString str=atompair.value.c_str();
            glPushMatrix();
            glTranslatef ( middle.x(),middle.y()-0.25,middle.z() );
            GLVisorBase::RenderText ( 0.0,0.0,0.0,str );
            glPopMatrix();


            glDisable ( GL_LINE_STIPPLE );
            glEnable ( GL_LIGHTING );
        }


    }

    if (m_banglechange)
    {

        Molecule::triad atompair = m_angles [m_pair];

        static float color[4];
        color[0]=0.9f;
        color[1]=0.3f;
        color[2]=0.3f;
        color[3]=0.5f;

        glMaterialfv ( GL_FRONT_AND_BACK, GL_AMBIENT, color );
        glMaterialfv ( GL_FRONT_AND_BACK, GL_DIFFUSE, color );

        glPushMatrix();
        glTranslatef ( molecule.CurrentFrame().XYZ()[atompair.i].x(),
                molecule.CurrentFrame().XYZ()[atompair.i].y(),
                molecule.CurrentFrame().XYZ()[atompair.i].z()  );
        gluSphere(quadric, 0.25, 8, 8);
        glPopMatrix();
        glPushMatrix();
        glTranslatef( molecule.CurrentFrame().XYZ()[atompair.j].x(),
                molecule.CurrentFrame().XYZ()[atompair.j].y(),
                molecule.CurrentFrame().XYZ()[atompair.j].z() );
        gluSphere(quadric, 0.25, 8, 8);
        glPopMatrix();
        glPushMatrix();
        glTranslatef( molecule.CurrentFrame().XYZ()[atompair.k].x(),
                molecule.CurrentFrame().XYZ()[atompair.k].y(),
                molecule.CurrentFrame().XYZ()[atompair.k].z() );
        gluSphere(quadric, 0.25, 8, 8);
        glPopMatrix();
    }

    if (m_bdihedralchange)
    {
        Molecule::quad atompair = m_dihedrals [m_pair];

        static float color[4];
        color[0]=0.9f;
        color[1]=0.3f;
        color[2]=0.3f;
        color[3]=0.5f;

        glMaterialfv ( GL_FRONT_AND_BACK, GL_AMBIENT, color );
        glMaterialfv ( GL_FRONT_AND_BACK, GL_DIFFUSE, color );

        glPushMatrix();
        glTranslatef ( molecule.CurrentFrame().XYZ()[atompair.i].x(),
                molecule.CurrentFrame().XYZ()[atompair.i].y(),
                molecule.CurrentFrame().XYZ()[atompair.i].z()  );
        gluSphere(quadric, 0.25, 8, 8);
        glPopMatrix();
        glPushMatrix();
        glTranslatef( molecule.CurrentFrame().XYZ()[atompair.j].x(),
                molecule.CurrentFrame().XYZ()[atompair.j].y(),
                molecule.CurrentFrame().XYZ()[atompair.j].z() );
        gluSphere(quadric, 0.25, 8, 8);
        glPopMatrix();
        glPushMatrix();
        glTranslatef( molecule.CurrentFrame().XYZ()[atompair.k].x(),
                molecule.CurrentFrame().XYZ()[atompair.k].y(),
                molecule.CurrentFrame().XYZ()[atompair.k].z() );
        gluSphere(quadric, 0.25, 8, 8);
        glPopMatrix();
        glPushMatrix();
        glTranslatef( molecule.CurrentFrame().XYZ()[atompair.l].x(),
                molecule.CurrentFrame().XYZ()[atompair.l].y(),
                molecule.CurrentFrame().XYZ()[atompair.l].z() );
        gluSphere(quadric, 0.25, 8, 8);
        glPopMatrix();

    }

    const int spheres=16;
    if ( m_bshowdipole )
    {
        if ( molecule.CurrentFrame().GetDipole().Norm() > 0.1 )
        {
            static float color[4];
            color[0]=0.4;
            color[1]=0.4;
            color[2]=0.9;
            color[3]=0.5;
            glMaterialfv ( GL_FRONT_AND_BACK, GL_AMBIENT, color );
            glMaterialfv ( GL_FRONT_AND_BACK, GL_DIFFUSE, color );
            glPushMatrix();
            {
                Coordinate zaxis ( 0,0,1 );
                Coordinate normal=zaxis^molecule.CurrentFrame().GetDipole();
                float angle=Coordinate::Angle( zaxis,molecule.CurrentFrame().GetDipole() );
                glRotatef ( angle*180/M_PI,normal.x(),normal.y(),normal.z() );
                float height=molecule.CurrentFrame().GetDipole().Norm()*0.3;
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

                    QColor col ( QColor ( 255,255,255 ) );
                    qglColor ( col );

                    //write two spaces to leave some separation between
                    //number and atom
                    QString dipole;
                    dipole.sprintf ( "  %.3f debyes",molecule.CurrentFrame().GetDipole().Norm());


                    RenderText(0.0,0.0,0.0,dipole);//RenderText ( dipole,0,0, 0 );
                    glEnable ( GL_LIGHTING );
                }
                glPopMatrix();
            }
            glPopMatrix();
        }
    }

    if ( m_bshowcell )
    {
        //                if ( !m_molecule.Cell().empty() )
        //                {
        //                  if ( mode == GL_RENDER )
        //                  {
        //                    glPushMatrix();
        //                    Coordinate c=m_molecule.Cell().at ( 0 ) +m_molecule.Cell().at ( 1 ) +m_molecule.Cell().at ( 2 );
        //                    c=-0.5*c;
        //                    glTranslatef ( c.x,c.y,c.z );
        //                    glDisable ( GL_LIGHTING );
        //                    std::cout << "show cell " << c << std::endl;
        //                    std::vector<Coordinate> cube=CalculateCellCube();
        //                    glColor4f ( 0.7,0.7,0.95,0.15 );
        //                    glBegin ( GL_QUADS );
        //                    std::cout << cube[0] << "," << cube[1] << "," << cube[2] << std::endl;
        //                    glVertexC ( cube[0] );
        //                    glVertexC ( cube[1] );
        //                    glVertexC ( cube[2] );
        //                    glVertexC ( cube[3] );

        //                    glVertexC ( cube[7] );
        //                    glVertexC ( cube[6] );
        //                    glVertexC ( cube[5] );
        //                    glVertexC ( cube[4] );

        //                    glVertexC ( cube[1] );
        //                    glVertexC ( cube[2] );
        //                    glVertexC ( cube[6] );
        //                    glVertexC ( cube[5] );

        //                    glVertexC ( cube[0] );
        //                    glVertexC ( cube[3] );
        //                    glVertexC ( cube[7] );
        //                    glVertexC ( cube[4] );

        //                    glVertexC ( cube[2] );
        //                    glVertexC ( cube[3] );
        //                    glVertexC ( cube[7] );
        //                    glVertexC ( cube[6] );

        //                    glVertexC ( cube[0] );
        //                    glVertexC ( cube[1] );
        //                    glVertexC ( cube[5] );
        //                    glVertexC ( cube[4] );


        //                    glEnd();
        //                    glEnable ( GL_LIGHTING );
        //                    glPopMatrix();
        //                  }
        //                }
    }

    if (m_bshowdensity)
    {
        const std::vector<RenderDensity>& positiverenderdensityvector = m_world->CurrentMolecule()->CurrentFrame().PositiveDensity();

        glPushMatrix();

        static float vcolor[4];
        vcolor[0]=0.0;
        vcolor[1]=0.1;
        vcolor[2]=0.9;
        vcolor[3]=m_transparence;
        glMaterialfv ( GL_FRONT_AND_BACK, GL_AMBIENT, vcolor );
        glMaterialfv ( GL_FRONT_AND_BACK, GL_DIFFUSE, vcolor );
        //glPushMatrix();
        glBegin(GL_TRIANGLES);
        for (size_t i=0; i<positiverenderdensityvector.size(); i++)
        {
            const RenderDensity& prv=positiverenderdensityvector[i];
            GLint flagindex=prv.FlagIndex();

            //Draw the triangles that were found.  There can be up to five per cube
            for (GLint triangle = 0; triangle < 5; triangle++)
            {
                if (triangleConnectionTable[flagindex][3*triangle] < 0) break;

                for (GLint corner = 0; corner < 3; corner++)
                {
                    GLint vertex = triangleConnectionTable[flagindex][3*triangle+corner];

                    const RenderDensity::GLcoordinate& edgeVertex = prv.EdgeVertex()[vertex];
                    const RenderDensity::GLcoordinate& edgeNorm = prv.EdgeNorm()[vertex];
                    {
                        glNormal3f(edgeNorm.X,edgeNorm.Y,edgeNorm.Z);
                        glVertex3f(edgeVertex.X,edgeVertex.Y,edgeVertex.Z);
                    }
                }
            }

        }
        glEnd();
        //glPopMatrix();
        glPopMatrix();

        const std::vector<RenderDensity>& negativerenderdensityvector = m_world->CurrentMolecule()->CurrentFrame().NegativeDensity();

        glPushMatrix();

        vcolor[0]=0.9;
        vcolor[1]=0.1;
        vcolor[2]=0.0;
        vcolor[3]=m_transparence;
        glMaterialfv ( GL_FRONT_AND_BACK, GL_AMBIENT, vcolor );
        glMaterialfv ( GL_FRONT_AND_BACK, GL_DIFFUSE, vcolor );
        //glPushMatrix();
        glBegin(GL_TRIANGLES);
        for (size_t i=0; i<negativerenderdensityvector.size(); i++)
        {
            const RenderDensity& nrv=negativerenderdensityvector[i];
            GLint flagindex=nrv.FlagIndex();

            //Draw the triangles that were found.  There can be up to five per cube
            for (GLint triangle = 0; triangle < 5; triangle++)
            {
                if (triangleConnectionTable[flagindex][3*triangle] < 0) break;

                for (GLint corner = 0; corner < 3; corner++)
                {
                    GLint vertex = triangleConnectionTable[flagindex][3*triangle+corner];

                    const RenderDensity::GLcoordinate& edgeVertex = nrv.EdgeVertex()[vertex];
                    const RenderDensity::GLcoordinate& edgeNorm = nrv.EdgeNorm()[vertex];
                    {
                        glNormal3f(edgeNorm.X,edgeNorm.Y,edgeNorm.Z);
                        glVertex3f(edgeVertex.X,edgeVertex.Y,edgeVertex.Z);
                    }
                }
            }

        }
        glEnd();
        //glPopMatrix();
        glPopMatrix();
    }


    if ( m_bshowaxis )
    {
        glPushMatrix();
        {
            static float color[4];
            color[0]=1.0;
            color[1]=0.0;
            color[2]=0.0;
            color[3]=1.0;
            glMaterialfv ( GL_FRONT_AND_BACK, GL_AMBIENT, color );
            glMaterialfv ( GL_FRONT_AND_BACK, GL_DIFFUSE, color );
            glPushMatrix();
            {
                Coordinate zaxis ( 0,0,1 );
                Coordinate axis (1,0,0);
                Coordinate normal=zaxis^axis;
                float angle=Coordinate::Angle( zaxis,axis );
                glRotatef ( angle*180/M_PI,normal.x(),normal.y(),normal.z() );
                gluCylinder ( quadric,0.05,0.05,1.0,spheres,spheres );
                QColor col ( QColor ( 255,0,0 ) );
                qglColor ( col );
                QString dipole;
                dipole.sprintf ( "X");
                RenderText(0.0,0.0,1.3,dipole);
                glPushMatrix();
                {
                    glTranslatef ( 0,0,1.0 );
                    gluCylinder ( quadric,0.09,0.001,0.2,spheres,spheres );
                }
                glPopMatrix();
            }
            glPopMatrix();

            color[0]=0.0;
            color[1]=1.0;
            color[2]=0.0;
            color[3]=1.0;
            glMaterialfv ( GL_FRONT_AND_BACK, GL_AMBIENT, color );
            glMaterialfv ( GL_FRONT_AND_BACK, GL_DIFFUSE, color );
            glPushMatrix();
            {
                Coordinate zaxis ( 0,0,1 );
                Coordinate axis (0,1,0);
                Coordinate normal=zaxis^axis;
                float angle=Coordinate::Angle( zaxis,axis );
                glRotatef ( angle*180/M_PI,normal.x(),normal.y(),normal.z() );
                gluCylinder ( quadric,0.05,0.05,1.0,spheres,spheres );
                QColor col ( QColor ( 0,255,0 ) );
                qglColor ( col );
                QString dipole;
                dipole.sprintf ( "Y");
                RenderText(0.0,0.0,1.3,dipole);

                glPushMatrix();
                {
                    glTranslatef ( 0,0,1.0 );
                    gluCylinder ( quadric,0.09,0.001,0.2,spheres,spheres );
                }
                glPopMatrix();
            }
            glPopMatrix();

            color[0]=0.0;
            color[1]=0.0;
            color[2]=1.0;
            color[3]=1.0;
            glMaterialfv ( GL_FRONT_AND_BACK, GL_AMBIENT, color );
            glMaterialfv ( GL_FRONT_AND_BACK, GL_DIFFUSE, color );
            glPushMatrix();
            {
                Coordinate zaxis ( 0,0,1 );
                Coordinate normal=zaxis^zaxis;
                float angle=Coordinate::Angle( zaxis,zaxis );
                glRotatef ( angle*180/M_PI,normal.x(),normal.y(),normal.z() );
                gluCylinder ( quadric,0.05,0.05,1.0,spheres,spheres );
                QColor col ( QColor ( 0,0,255 ) );
                qglColor ( col );
                QString dipole;
                dipole.sprintf ( "Z");
                RenderText(0.0,0.0,1.3,dipole);

                glPushMatrix();
                {
                    glTranslatef ( 0,0,1.0 );
                    gluCylinder ( quadric,0.09,0.001,0.2,spheres,spheres );
                }
                glPopMatrix();
            }
            glPopMatrix();
        }
        glPopMatrix();
    }

    gluDeleteQuadric ( quadric );
    glPopMatrix();

    return;

}


/** Put trackball rotation center on centroid of the visible part*/
void GLVisor::CenterVisiblePart()
{
    if ( !m_world->CurrentMolecule() ) return;
    if ( m_world->CurrentMolecule()->Residues().empty() ) return;
    const std::vector<Frame>& frames=m_world->CurrentMolecule()->Frames();
    Coordinate current;
    for(size_t i=0;i<frames.size();++i)
    {

        Coordinate c ( 0,0,0 );
        size_t natoms=m_world->CurrentMolecule()->Atoms().size();
        size_t visibleatoms=0;
        for ( size_t j=0;j<natoms;++j )
        {
            if ( m_world->CurrentMolecule()->Atoms() [j].Residue()->Visible() )
            {
                visibleatoms++;
                c+= frames[i].XYZ() [j];
            }
        }
        c/=visibleatoms;
        Handlers() [ m_world->CurrentMoleculeIndex() ].SetRotationCenter ( c,i );
        if ( i == m_world->CurrentMolecule()->CurrentFrameIndex() )
            current=c;
    }

    SetCamera ( current );
    SetupProjection();
    update();///

}

/** \brief Center molecule

if \a b is true (default) the whole molecule is centered otherwise
only visible part would be centered
*/

void GLVisor::CenterMolecule ( bool b/*=true*/, bool bblock/*=true*/ )
{
    if ( b )
    {
        if ( !m_world->Molecules().empty() )
        {
            std::vector<Frame>& frames=m_world->CurrentMolecule()->Frames();
            m_world->Visor()->Handlers()[m_world->CurrentMoleculeIndex()].SetBlocked(bblock);
            for ( size_t i=0;i<frames.size();++i )
            {

                Coordinate centroid=m_world->CurrentMolecule()->Frames()[i].Centroid();

                Handlers() [m_world->CurrentMoleculeIndex()].SetRotationCenter ( centroid,i );
            }
            SetCamera ( m_world->CurrentMolecule()->CurrentFrame().Centroid() );
            SetupProjection();
            update();
        }
    }
    else {
        CenterVisiblePart();
    }
}

/** Center current molecule*/
void GLVisor::Center()
{
    SetCamera ( Handlers()[m_world->CurrentMoleculeIndex()].RotationCenter(m_world->CurrentMolecule()->CurrentFrameIndex() ) );
    SetupProjection();
    update();
}

/** manage mouse wheel events*/
void GLVisor::wheelEvent ( QWheelEvent* e )
{

    if (!m_world->CurrentMolecule() ) return;

    if ( m_selmode == ROTATEBOND )
    {
        if ( m_selmode == ROTATEBOND && m_selatoms.size() ==2 )
        {
            m_world->CurrentMolecule()->RotateBond ( m_selatoms[0],m_selatoms[1],-e->delta() );
            update();
        }
    }
    else
    {
        if ( e->delta() <0 )
        {
            Handlers() [m_world->CurrentMoleculeIndex() ].Scale() *=1.3;
        }
        if ( e->delta() >0 )
        {
            Handlers() [m_world->CurrentMoleculeIndex() ].Scale() /=1.3;
        }
    }

    update();
}

void GLVisor::ResetTrackball()
{
    for(std::vector<MoleculeHandler>::iterator it=Handlers().begin();it!=Handlers().end();++it)
    {
        it->Reset();
    }
}

QStringList GLVisor::GetDistances (size_t frame)
{
    std::vector<Molecule::pair>::iterator it;
    const Frame& fr=m_world->CurrentMolecule()->Frames().at(frame);
    QStringList list;
    for (it=m_distances.begin();it!=m_distances.end();it++)
    {
        QString str;
        str.sprintf ( "(%s%d,%s%d)  %.3f",
                      m_world->CurrentMolecule()->Atoms() [it->i].Symbol().c_str(),
                it->i+1,
                m_world->CurrentMolecule()->Atoms() [it->j].Symbol().c_str(),
                it->j+1,
                fr.Distance ( it->i,it->j ) );
        it->value = str.toStdString();
        list << str;
    }
    return list;
}

QStringList GLVisor::GetAngles (size_t frame)
{
    std::vector<Molecule::triad>::iterator it;
    const Frame& fr=m_world->CurrentMolecule()->Frames().at(frame);
    QStringList list;
    for (it=m_angles.begin();it!=m_angles.end();it++)
    {
        QString str;
        str.sprintf ( "(%s%d,%s%d,%s%d)  %.3f",
                      m_world->CurrentMolecule()->Atoms() [it->i].Symbol().c_str(),
                it->i+1,
                m_world->CurrentMolecule()->Atoms() [it->j].Symbol().c_str(),
                it->j+1,
                m_world->CurrentMolecule()->Atoms() [it->k].Symbol().c_str(),
                it->k+1,
                fr.Angle ( it->i,it->j,it->k ) );
        it->value = str.toStdString();
        list << str;
    }
    return list;
}

QStringList GLVisor::GetDihedrals (size_t frame)
{
    std::vector<Molecule::quad>::iterator it;
    const Frame& fr=m_world->CurrentMolecule()->Frames().at(frame);
    QStringList list;
    for (it=m_dihedrals.begin();it!=m_dihedrals.end();it++)
    {
        QString str;
        str.sprintf ( "(%s%d,%s%d,%s%d,%s%d)  %.3f",
                      m_world->CurrentMolecule()->Atoms() [it->i].Symbol().c_str(),
                it->i+1,
                m_world->CurrentMolecule()->Atoms() [it->j].Symbol().c_str(),
                it->j+1,
                m_world->CurrentMolecule()->Atoms() [it->k].Symbol().c_str(),
                it->k+1,
                m_world->CurrentMolecule()->Atoms() [it->l].Symbol().c_str(),
                it->l+1,
                fr.Dihedral ( it->i,it->j,it->k,it->l )*180/M_PI );
        it->value = str.toStdString();
        list << str;
    }
    return list;
}

void GLVisor::OnClearMeasures()
{
    OnResetSelection();

    m_bshowdistances = false;
    m_bdistancechange = false;
    m_banglechange = false;
    m_bdihedralchange = false;

    m_distances.clear();
    m_angles.clear();
    m_dihedrals.clear();

    update();
}

void GLVisor::OnDistanceChange(int i)
{
    OnResetSelection();

    m_pair = i;

    m_bdistancechange = true;
    m_banglechange = false;
    m_bdihedralchange = false;

    update();
}

void GLVisor::OnAngleChange(int i)
{

    OnResetSelection();

    m_pair = i;

    m_bdistancechange = false;
    m_banglechange = true;
    m_bdihedralchange = false;

    update();

}

void GLVisor::OnDihedralChange(int i)
{
    OnResetSelection();

    m_pair = i;

    m_bdistancechange = false;
    m_banglechange = false;
    m_bdihedralchange = true;

    update();

}

void GLVisor::OnShowDistances(bool b)
{
    m_bshowdistances = b;

    update();
}

void GLVisor::OnShowDipole(bool b)
{
    m_bshowdipole=b;

    update();
}

void GLVisor::OnShowCell(bool b)
{
    m_bshowcell=b;

    update();
}

void GLVisor::OnShowAxis(bool b)
{
    m_bshowaxis=b;

    update();
}

void GLVisor::OnShowDensity(bool b)
{
    m_bshowdensity=b;

    update();
}

void GLVisor::OnTransparenceChange(float f)
{
    m_transparence=f;

    update();
}

void GLVisor::OnEnantiomerize()
{
    std::vector<Coordinate>::iterator it;
    for ( it=m_world->CurrentMolecule()->CurrentFrame().XYZ().begin();it!=m_world->CurrentMolecule()->CurrentFrame().XYZ().end();++it )
    {
        (*it) *=-1;
    }
    size_t i,j;
    for ( i=0;i<m_world->CurrentMolecule()->CurrentFrame().GetForces().size();++i )
    {
        m_world->CurrentMolecule()->CurrentFrame().GetForces() ( i ) *=-1;
    }
    for ( i=0;i<m_world->CurrentMolecule()->CurrentFrame().GetHessian().NRows();++i )
        for ( j=0;j<m_world->CurrentMolecule()->CurrentFrame().GetHessian().NColumns();++j )
            m_world->CurrentMolecule()->CurrentFrame().GetHessian() ( i,j ) *=-1;

    update();
}
