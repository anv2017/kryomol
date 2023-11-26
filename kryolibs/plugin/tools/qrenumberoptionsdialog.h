/*****************************************************************************************
                            qrenumberoptionsdialog.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QRENUMBEROPTIONSDIALOG_H
#define QRENUMBEROPTIONSDIALOG_H

#include <qdialog.h>

class QRenumberOptionsDialog : public QDialog
{
    Q_OBJECT

    public:
            QRenumberOptionsDialog ( QWidget* parent =0 );
            ~QRenumberOptionsDialog();
            enum RenumberMode { AllAtoms, HeavyAtoms };
            RenumberMode Mode() const { return m_rmode; }
    private slots:
            void OnHeavyAtoms();
            void OnAllAtoms();
    private:
            RenumberMode m_rmode;

};

#endif
