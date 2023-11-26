/*****************************************************************************************
                            qmolecularlistcontrol.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <fstream>
#include <sstream>

#include <QLineEdit>
#include <QValidator>
#include <QVBoxLayout>
#include <QFileDialog>
#include <QMessageBox>
#include <QClipboard>
#include <QButtonGroup>

#include "world.h"
#include "thermo.h"
#include "molecule.h"
#include "stringtools.h"

#include "qmolecularlistcontrol.h"
#include "qmoleculartreewidget.h"

using namespace kryomol;

QMolecularListControl::QMolecularListControl ( QWidget* parent ) : QWidget ( parent )
{  
    setupUi ( this );
    m_ensamblegroup = new QButtonGroup ( this );
    m_ensamblegroup->addButton ( _nveRadioButton,0 );
    m_ensamblegroup->addButton ( _nvtRadioButton,1 );
    QVBoxLayout* vbox= new QVBoxLayout ( _listFrame );
    m_listview = new QMolecularTreeWidget ( _listFrame );
    vbox->addWidget ( m_listview );
    QStringList initiallabel;
    initiallabel << "Conformers";
    m_listview->setHeaderLabels ( initiallabel );


    connect ( m_ensamblegroup,SIGNAL ( buttonClicked ( int ) ),this,SLOT ( OnChangeEnsamble ( int ) ) );

    QValidator* tval= new QDoubleValidator ( 0.001,30000,3,_TLineEdit );
    _TLineEdit->setValidator ( tval );
    connect ( _TLineEdit,SIGNAL ( returnPressed() ),this,SLOT ( OnChangeTemperature() ) );
    connect(_calculateButton,SIGNAL(clicked()),this,SLOT(OnCalculateButton()));
    connect(_readEnergiesButton,SIGNAL(clicked()),this,SLOT(OnReadEnergiesButton()));
    connect(_readPopulationsButton,SIGNAL(clicked()),this,SLOT(OnReadPopulationsButton()));
    connect(_pasteEnergiesButton,SIGNAL(clicked()),this,SLOT(OnPasteEnergiesButton()));
    connect(_pastePopulationsButton,SIGNAL(clicked()),this,SLOT(OnPastePopulationsButton()));


}


QMolecularListControl::~QMolecularListControl()
{}


void QMolecularListControl::Init()
{
    if ( !m_world ) return;
    QString str;
    str.sprintf ( "%.2f",m_world->Temperature() );
    _TLineEdit->setText ( str );
    m_listview->clear();
    m_listview->Init();
}

void QMolecularListControl::SetWorld ( World* world )
{
    m_world=world;
    m_listview->SetWorld ( world );
    QString ts;
    ts.sprintf ( "%.3f",m_world->Temperature() );
    _TLineEdit->setText ( ts );
    m_ensamblegroup->button ( static_cast<int> ( m_world->Ensamble() ) )->setChecked ( true );
}

void QMolecularListControl::OnChangeEnsamble ( int id )
{
    m_world->SetEnsamble ( static_cast<kryomol::ensamble> ( id ) );
}

void QMolecularListControl::OnChangeTemperature()
{
    m_world->SetTemperature ( _TLineEdit->text().toDouble() );

}

void QMolecularListControl::OnCalculateButton()
{
    try
    {
        m_world->SetTemperature ( _TLineEdit->text().toDouble() );
        m_world->SetPopulations(*m_world->CurrentMolecule());
    }
    catch(std::exception& e)
    {
        QMessageBox::critical(this,"",QString(e.what()));
    }
    m_listview->Init();
}

void QMolecularListControl::ReadPopulationsFromStream(std::istream& stream)
{


    std::vector<double> pops;
    std::string line,wholeline;
    while(std::getline(stream,line))
    {
        line+="\t";
        wholeline+=line;
    }

    StringTokenizer tok(wholeline," ,\t\r");
    for(StringTokenizer::iterator it=tok.begin();it!=tok.end();++it)
    {
        float pop=kryomol::atof(*it);
        if (  pop < 0 ) throw kryomol::Exception("Populations should be real positive numbers");
        pops.push_back(kryomol::atof(*it));
    }


    if ( pops.size() != m_world->CurrentMolecule()->Frames().size() )
        throw kryomol::Exception("The number of populations is not the same as the number of conformers in current molecule");

    //normalize the populations
    double sum=0;
    for(size_t i=0;i<pops.size();++i)
        sum+=pops[i];
    for(size_t i=0;i<pops.size();++i)
        pops[i]/=sum;
    m_world->CurrentMolecule()->Populations()=pops;
}

void QMolecularListControl::ReadEnergiesFromStream(std::istream& stream)
{
    std::vector<double> energies;

    std::string line,wholeline;
    while(std::getline(stream,line))
    {
        line+="\t";
        wholeline+=line;
    }

    StringTokenizer tok(wholeline," ,\t\r");
    if ( !tok.empty() )
    {
        for(StringTokenizer::iterator it=tok.begin();it!=tok.end();++it)
        {
            energies.push_back(kryomol::atof(*it));
        }
    }

    if ( energies.size() != m_world->CurrentMolecule()->Frames().size() )
        throw kryomol::Exception("The number of energies read is not the same as the number of conformers in current molecule");

    for(size_t i=0;i<energies.size();++i)
    {
        m_world->CurrentMolecule()->Frames().at(i).PotentialEnergy()=energies[i];
    }

}

void QMolecularListControl::OnReadEnergiesButton()
{
    QString filename=QFileDialog::getOpenFileName(this,tr("Read Energies"));
    if ( !filename.isEmpty() )
    {
        std::ifstream fstream(filename.toLatin1());
        this->ReadEnergiesFromStream(fstream);
    }
    m_listview->Init();
}

void QMolecularListControl::OnReadPopulationsButton()
{
    QString filename=QFileDialog::getOpenFileName(this,tr("Read Populations"));
    if ( !filename.isEmpty() )
    {
        std::ifstream fstream(filename.toLatin1());
        this->ReadPopulationsFromStream(fstream);
    }
    m_listview->Init();
}

void QMolecularListControl::OnPasteEnergiesButton()
{

    QClipboard *clipboard = QApplication::clipboard();
    if ( clipboard )
    {
        QString text=clipboard->text();
        if ( !text.isEmpty() )
        {
            std::stringstream stream;
            stream << text.toStdString();
            this->ReadEnergiesFromStream(stream);
        }
    }

    m_listview->Init();

}

void QMolecularListControl::OnPastePopulationsButton()
{
    QClipboard *clipboard = QApplication::clipboard();
    if ( clipboard )
    {
        QString text=clipboard->text();
        if ( !text.isEmpty() )
        {
            std::stringstream stream;
            stream << text.toStdString();
            this->ReadPopulationsFromStream(stream);
        }
    }
    m_listview->Init();

}


