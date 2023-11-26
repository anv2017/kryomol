/*****************************************************************************************
                            qnumlistviewitem.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QNUMLISTVIEWITEM_H
#define QNUMLISTVIEWITEM_H

#include <QTreeWidgetItem>
#include "export.h"

namespace kryomol
{
  class KRYOMOL_API QNumTreeWidgetItem : public QTreeWidgetItem
  {
    public:
      QNumTreeWidgetItem ( int type = Type );
      QNumTreeWidgetItem ( const QStringList & strings, int type = Type );
      QNumTreeWidgetItem ( QTreeWidget * parent, int type = Type );
      QNumTreeWidgetItem ( QTreeWidget * parent, const QStringList & strings, int type = Type );
      QNumTreeWidgetItem ( QTreeWidget * parent, QTreeWidgetItem * preceding, int type = Type );
      QNumTreeWidgetItem ( QTreeWidgetItem * parent, int type = Type );
      QNumTreeWidgetItem ( QTreeWidgetItem * parent, const QStringList & strings, int type = Type );
      QNumTreeWidgetItem ( QTreeWidgetItem * parent, QTreeWidgetItem * preceding, int type = Type );   
      bool operator< ( const QTreeWidgetItem &other ) const;


  };
}
#endif
