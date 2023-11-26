/*****************************************************************************************
                            qpdbcontrol.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/


#include <QTreeWidget>
#include <QLayout>
#include <QMenu>

#include <qnumlistviewitem.h>

#include "qpdbcontrol.h"
#include "pdbtools.h"
#include "molecule.h"
#include "world.h"
#include "glvisor.h"

using namespace kryomol;

QPDBControl::QPDBControl ( QWidget* parent ) : QWidget ( parent )
{
	QVBoxLayout* lay = new QVBoxLayout ( this );
	m_tree= new QPDBTreeWidget ( this );
	lay->addWidget ( m_tree );

}

void QPDBControl::Init()
{
  m_tree->clear();
	m_tree->Init();

}

QPDBTreeWidget::QPDBTreeWidget ( QWidget* parent ) : QTreeWidget ( parent )
{
	m_menu = new QMenu ( this );
	QAction* showAction = new QAction ( "Show",this );
	QAction* hideAction = new QAction ( "Hide",this );
	connect ( showAction,SIGNAL ( activated() ),this,SLOT ( OnShowAction() ) );
	connect ( hideAction,SIGNAL ( activated() ),this,SLOT ( OnHideAction() ) );
        QStringList label;
        label << "No PDB Info";
        setHeaderLabels(label);

}
void QPDBTreeWidget::Init()
{
	const std::vector<PDBResidue*>& vres=m_world->CurrentMolecule()->Residues();
	size_t residues=vres.size();

	this->setColumnCount ( 2 );
	QStringList labels;
	labels << "Residue" << "Number";
	setHeaderLabels ( labels );
	QTreeWidgetItem* root= new QTreeWidgetItem ( this );
	root->setFlags ( Qt::ItemIsSelectable |   Qt::ItemIsUserCheckable | Qt::ItemIsEnabled );

	root->setText ( 0,"Molecule" );
	this->addTopLevelItem ( root );

	for ( size_t i=0;i<residues;i++ )
	{
		QString res ( vres[i]->Name().c_str() );
		QString nres(vres[i]->Index().c_str());
		QStringList sl;
                sl << res << nres;
		QNumTreeWidgetItem* rit= new QNumTreeWidgetItem ( root,sl );
		rit->setFlags ( Qt::ItemIsSelectable |   Qt::ItemIsUserCheckable | Qt::ItemIsEnabled );
		if ( vres[i]->Visible() )
			rit->setCheckState ( 0,Qt::Checked );
		else
			rit->setCheckState ( 0,Qt::Unchecked );

	}
	root->setCheckState ( 0,Qt::Checked );
	setSelectionMode ( QAbstractItemView::ExtendedSelection );
	connect ( this,SIGNAL ( itemChanged ( QTreeWidgetItem*,int ) ),this,SLOT ( OnItemChanged ( QTreeWidgetItem* ,int ) ) );
}

void QPDBTreeWidget::OnShowAction()
{
	SetVisibility ( true );
}

void QPDBTreeWidget::OnHideAction()
{
	SetVisibility ( false );
}

void QPDBTreeWidget::SetVisibility ( bool b )
{
	QTreeWidgetItem* root=this->topLevelItem ( 0 );
	int nchilds=root->childCount();
	for ( int i=0;i<nchilds;++i )
	{
		m_world->CurrentMolecule()->Residues().at ( i )->SetVisible ( b );
	}
	m_world->Visor()->update();
}

void QPDBTreeWidget::OnItemChanged ( QTreeWidgetItem* item,int column )
{
	if ( column != 0 ) return;
	if ( item == this->topLevelItem ( 0 ) )
	{
		bool checked=item->checkState ( 0 ) == Qt::Checked;
		SetVisibility ( checked );
		disconnect ( this,SIGNAL ( itemChanged ( QTreeWidgetItem*,int ) ),this,SLOT ( OnItemChanged ( QTreeWidgetItem* ,int ) ) );
		int nchilds=item->childCount();
		for ( int i=0;i<nchilds;++i )
		{
			if ( checked )
				item->child ( i )->setCheckState ( 0,Qt::Checked );
			else
				item->child ( i )->setCheckState ( 0,Qt::Unchecked );
		}
		connect ( this,SIGNAL ( itemChanged ( QTreeWidgetItem*,int ) ),this,SLOT ( OnItemChanged ( QTreeWidgetItem* ,int ) ) );
		return;
        }
	std::string name=item->text ( 0 ).toStdString();
	std::string index=item->text ( 1 ).toStdString();
	std::vector<PDBResidue*>::iterator rt;

	for ( rt=m_world->CurrentMolecule()->Residues().begin();rt!=m_world->CurrentMolecule()->Residues().end();++rt )
	{
		if ( ( *rt )->Name() == name && ( *rt )->Index() == index ) break;
	}
	if ( rt == m_world->CurrentMolecule()->Residues().end() )
	{
		std::cerr << "Cannot determine residue" << std::endl;
		return;
	}
	if ( item->checkState ( 0 ) == Qt::Checked )
	{

		( *rt )->SetVisible ( true );
	}
	else
        {
		( *rt )->SetVisible ( false );
	}

	m_world->Visor()->update();
}
