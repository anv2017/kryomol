/*****************************************************************************************
                            qjoboptwidget.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QJOBOPTWIDGET_H
#define QJOBOPTWIDGET_H

#include <QWidget>

#include "qjobwidget.h"
#include "qconvwidget.h"

namespace kryomol
{
class World;

};

class QJobOptWidget : public QJobWidget
{

    Q_OBJECT

public:
    QJobOptWidget(kryomol::World* world, QWidget* parent = 0);
    ~QJobOptWidget();

private:
    void Init();

private:
    QConvWidget* m_convwidget;

};


#endif // QJOBOPTWIDGET_H
