/*****************************************************************************************
                            qmoleculartreewidget.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QMOLECULARTREEWIDGET_H
#define QMOLECULARTREEWIDGET_H

#include <QTreeWidget>
#include <QMouseEvent>

/**
A specialized list view for molecular data
*/

namespace kryomol
{
class GLVisor;
class World;
}

class QTreeWidgetItem;
class QMolecularTreeWidget : public QTreeWidget
{

Q_OBJECT
public:
  QMolecularTreeWidget(QWidget* parent=0);
  ~QMolecularTreeWidget();
  void SetWorld(kryomol::World* w);
  void Init();
signals:
  void onItem(QTreeWidgetItem* );
private slots:
  void OnItemChanged(QTreeWidgetItem* item);
  void OnFrame(size_t );

private: 
  kryomol::World* m_world;
  int m_current;

};

#endif
