#include "qbasejobwidget.h"
#include "nmrdevapp.h"
#include "world.h"
#include "atom.h"
#include "glvisor.h"

#include <QMainWindow>


QBaseJobWidget::QBaseJobWidget ( QWidget* parent) : QWidget(parent)
{
    setAttribute(Qt::WA_DeleteOnClose);
    Init();
}


QBaseJobWidget::QBaseJobWidget ( const char* file, QWidget* parent) : QWidget (parent) , m_file ( file )
{
  Init();
}

void QBaseJobWidget::Init()
{
    m_app=static_cast<nmrdev::NMRDevApplication*>qApp;

    m_world = new nmrdev::World(this);

    //m_world->Visor()->setSizePolicy ( QSizePolicy ( QSizePolicy::Expanding,QSizePolicy::Expanding ) );

    m_visor = new nmrdev::GLVisor(m_world);

    nmrdev::BuildPeriodicTable();

}

QBaseJobWidget::~QBaseJobWidget()
{

}


