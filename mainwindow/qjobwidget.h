
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
#include <QMainWindow>

namespace kryomol
{
class World;

};


class QDockWidget;
class QTabWidget;

class QOrbitalWidget;
class QJobWidget : public QMainWindow
{

    Q_OBJECT

public:
    QJobWidget(QWidget *parent = 0);
    ~QJobWidget();
    //virtual BuilldVisor();
    kryomol::World* World() {return m_world;}
    //void SetWorld(kryomol::World* world) {m_world = world;}
    QList<QDockWidget*> DockWidgets() const { return m_dockwidgets; }
    virtual void InitWidgets() {}
protected:
    void InitCommonWidgets();
    virtual void showEvent(QShowEvent *event);
    virtual void hideEvent(QHideEvent *event);

protected:
    kryomol::World* m_world;
    QList<QDockWidget*> m_dockwidgets;
    QTabWidget* m_tabwidget;

};

#endif // QJOBWIDGET_H
