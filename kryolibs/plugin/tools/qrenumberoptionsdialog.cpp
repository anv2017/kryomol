/*****************************************************************************************
                            qrenumberoptionsdialog.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <QPushButton>
#include <qlayout.h>
#include <qlabel.h>
#include "qrenumberoptionsdialog.h"
//Added by qt3to4:
#include <QGridLayout>

QRenumberOptionsDialog::QRenumberOptionsDialog ( QWidget* parent )
                : QDialog ( parent )
{
        QPushButton* allAtomsButton = new QPushButton ( this );
        QGridLayout* box = new QGridLayout ( this);
        QLabel* label= new QLabel(this);
        label->setText("Click atoms to renumber. Press Escape key to reset");
        allAtomsButton->setText ( "All atoms" );
        QPushButton* heavyAtomsButton = new QPushButton ( this );
        heavyAtomsButton->setText ( "Only heavy atoms" );
        QPushButton* cancelButton = new QPushButton ( this );
        cancelButton->setText ( "Cancel" );

        box->addWidget(label,0,0,0,2,Qt::AlignHCenter);
        box->addWidget ( allAtomsButton,1,0 );
        box->addWidget ( heavyAtomsButton,1,1 );
        box->addWidget ( cancelButton, 1,2 );
        connect ( allAtomsButton,SIGNAL ( clicked() ),this,SLOT ( OnAllAtoms() ) );
        connect ( heavyAtomsButton,SIGNAL ( clicked() ),this,SLOT ( OnHeavyAtoms() ) );
        connect ( cancelButton,SIGNAL ( clicked() ),this,SLOT ( reject() ) );

}


QRenumberOptionsDialog::~QRenumberOptionsDialog()
{}

void QRenumberOptionsDialog::OnAllAtoms()
{
        m_rmode=AllAtoms;
        accept();
}

void QRenumberOptionsDialog::OnHeavyAtoms()
{
        m_rmode= HeavyAtoms;
        accept();
}

