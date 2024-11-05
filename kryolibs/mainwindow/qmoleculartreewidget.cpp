/*****************************************************************************************
                            qmoleculartreewidget.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <QLabel>

#include "qmoleculartreewidget.h"

#include "world.h"
#include "qnumlistviewitem.h"
#include "thermo.h"
#include "molecule.h"

using namespace kryomol;

QMolecularTreeWidget::QMolecularTreeWidget ( QWidget* parent  )
    : QTreeWidget ( parent ) , m_world(NULL)
{
    connect ( this,SIGNAL ( currentItemChanged ( QTreeWidgetItem*, QTreeWidgetItem* ) ),this,SLOT ( OnItemChanged ( QTreeWidgetItem* ) ) );

}


QMolecularTreeWidget::~QMolecularTreeWidget()
{}


void QMolecularTreeWidget::Init()
{
    clear();
    QStringList labels;
    int ncolumns=2;
    labels << "#";
    labels << "color";
    const std::vector<double>& populations=m_world->CurrentMolecule()->Populations();

    if ( m_world->CurrentMolecule()->Frames().back().PotentialEnergy() )
    {
        labels << "Energies (a.u.)";
        ncolumns++;
    }

    if ( m_world->CurrentMolecule()->Populations().size() > 1 )
    {
        ncolumns++;
        labels << "Populations (%)";
    }
    setColumnCount(ncolumns);

    setHeaderLabels(labels);


    for ( size_t i=0;i<m_world->CurrentMolecule()->Frames().size();++i)
    {
        //m_visor->OnSelectPoint ( i-1,false );
        QString frame,energy, venergy, kenergy,spop;
        frame.sprintf ( "%d",static_cast<int> ( i+1 ) );

        double ve;
        const Frame& mf=m_world->CurrentMolecule()->Frames().at(i);
        QStringList sl;
        sl << frame;
        sl <<  QString();
        if ( mf.PotentialEnergy() )
        {
            ve=mf.PotentialEnergy().Value();;
            venergy.sprintf ( "%.5f",ve );
            sl << venergy;

        }

        if (!populations.empty() )
        {
            spop.sprintf("%.2f",populations.at(i)*100);
        }

        sl << spop;

        QNumTreeWidgetItem* mitem= new QNumTreeWidgetItem ( this, sl );
        float h,s,l;
        mf.GetColor(h,s,l);
        QFrame *f = new QFrame(this);
        QColor c=QColor::fromHslF(h,s,l);

        f->setStyleSheet(QString("QFrame { background-color: %1; border: none }").arg(c.name()));

        this->setItemWidget(mitem,1,f);
        addTopLevelItem(mitem);


    }

    setSortingEnabled(true);
    sortItems(0,Qt::AscendingOrder);
    for( int i=0;i<ncolumns;++i)
        this->resizeColumnToContents(i);

    OnFrame ( m_world->CurrentMolecule()->CurrentFrameIndex() );


}



void QMolecularTreeWidget::OnItemChanged ( QTreeWidgetItem* item )
{

    disconnect ( m_world,SIGNAL ( currentFrame ( size_t ) ),this,SLOT ( OnFrame ( size_t ) ) );
    if ( !item )
        return;
    if ( m_world )
    {
        int point =item->text ( 0 ).toInt();
        if ( point > 0 )
            m_world->SelectFrame( point -1 );
    }
    connect ( m_world,SIGNAL ( currentFrame ( size_t ) ),this,SLOT ( OnFrame ( size_t ) ) );
}


void QMolecularTreeWidget::SetWorld ( kryomol::World* w )
{
    m_world=w;
    connect ( m_world,SIGNAL ( currentFrame ( size_t ) ),this,SLOT ( OnFrame ( size_t ) ) );



}

void QMolecularTreeWidget::OnFrame ( size_t frame )
{

    int i=0;
    for(;i<topLevelItemCount();++i)
    {

        if ( ( atoi ( topLevelItem(i)->text ( 0 ).toStdString().c_str() )-1 ) == (int)frame ) break;

    }
    disconnect ( this,SIGNAL ( currentItemChanged ( QTreeWidgetItem*,QTreeWidgetItem* ) ),this,SLOT ( OnItemChanged( QTreeWidgetItem* ) ) );
    setCurrentItem ( topLevelItem( i ) );
    connect ( this,SIGNAL ( currentItemChanged ( QTreeWidgetItem*,QTreeWidgetItem* ) ),this,SLOT ( OnItemChanged ( QTreeWidgetItem* ) ) );
}


