/*****************************************************************************************
                            qpdbcontrol.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QPDBCONTROL_H
#define QPDBCONTROL_H
#include <QTreeWidget>

namespace kryomol
{
  class World;
}

class QMenu;
class QPDBTreeWidget : public QTreeWidget
{
  Q_OBJECT
  public:
  QPDBTreeWidget(QWidget* parent=0);
  void SetWorld(kryomol::World* w) { m_world=w; }
  void Init();
  private slots:
  void OnShowAction();
  void OnHideAction();
  void OnItemChanged(QTreeWidgetItem* item,int column);
  private:
  void SetVisibility(bool b);
  private:
  kryomol::World* m_world;
  QMenu* m_menu;
};

class QPDBControl : public QWidget
{
  public:
  QPDBControl(QWidget* parent=0);
  void SetWorld(kryomol::World* w) { m_tree->SetWorld(w); }
  void Init();
  private:
  QPDBTreeWidget* m_tree;

};

#endif
