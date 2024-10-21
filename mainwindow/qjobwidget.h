/*****************************************************************************************
                            qjobwidget.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QJOBWIDGET_H
#define QJOBWIDGET_H

#include <QWidget>
#include <QSplitter>


namespace kryomol
{
class World;

};


class QDockWidget;

class QJobWidget : public Q
{

    Q_OBJECT

public:
    QJobWidget(kryomol::World* world, QWidget *parent = 0);
    ~QJobWidget();
    virtual BuilldVisor();
    kryomol::World* GetWorld() {return m_world;}
    void SetWorld(kryomol::World* world) {m_world = world;}
    QList<QDockWidget*> DockWidgets() const { return m_dockwidgets; }

private:
    kryomol::World* m_world;
    QList<QWidget*> m_dockwidgets;

};

#endif // QJOBWIDGET_H
