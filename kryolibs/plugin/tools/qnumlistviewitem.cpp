/*****************************************************************************************
                            qnumlistviewitem.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include "qnumlistviewitem.h"

using namespace kryomol;

QNumTreeWidgetItem::QNumTreeWidgetItem ( int type ) : QTreeWidgetItem ( type )
{
}

QNumTreeWidgetItem::QNumTreeWidgetItem ( const QStringList & strings, int type ) : QTreeWidgetItem ( strings,type )
{
}

QNumTreeWidgetItem:: QNumTreeWidgetItem ( QTreeWidget * parent, int type ) : QTreeWidgetItem ( parent,type )
{
}

QNumTreeWidgetItem::QNumTreeWidgetItem ( QTreeWidget * parent, const QStringList & strings, int type ) : QTreeWidgetItem ( parent,strings,type )
{
}

QNumTreeWidgetItem::QNumTreeWidgetItem ( QTreeWidget * parent, QTreeWidgetItem * preceding, int type ) : QTreeWidgetItem ( parent,preceding,type )
{
}
QNumTreeWidgetItem::QNumTreeWidgetItem ( QTreeWidgetItem * parent, int type ) : QTreeWidgetItem ( parent,type )
{
}

QNumTreeWidgetItem::QNumTreeWidgetItem ( QTreeWidgetItem * parent, const QStringList & strings, int type ) : QTreeWidgetItem ( parent,strings,type )
{
}

QNumTreeWidgetItem::QNumTreeWidgetItem ( QTreeWidgetItem * parent, QTreeWidgetItem * preceding, int type ) : QTreeWidgetItem ( parent,preceding,type )
{
}

bool QNumTreeWidgetItem::operator< ( const QTreeWidgetItem &other ) const
{
  int column = treeWidget()->sortColumn();
  double a=this->text ( column ).toDouble();
  double b=other.text ( column ).toDouble();
  return ( a < b );
}
