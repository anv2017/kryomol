/*****************************************************************************************
                            qnicsgriddialog.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QNICSGRIDDIALOG_H
#define QNICSGRIDDIALOG_H

#include <QDialog>
#include <vector>
#include "coordinate.h"

namespace Ui {
class QNICSGridDialog;
}

class QNICSGridDialog : public QDialog
{
    Q_OBJECT
    
public:
    explicit QNICSGridDialog(std::pair<kryomol::Coordinate,kryomol::Coordinate>& box, QWidget *parent = 0);
    ~QNICSGridDialog();
    std::pair<kryomol::Coordinate,kryomol::Coordinate> Box() const;
    double GridResolution() const;

private:
    Ui::QNICSGridDialog *ui;
};

#endif // QNICSGRIDDIALOG_H
