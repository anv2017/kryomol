#ifndef CONFMANAGER_H
#define CONFMANAGER_H

#include <QWidget>

#include "kryomolcore_export.h"

namespace kryomol
{
class World;
}

class QTreeWidget;

class KRYOMOLCORE_EXPORT ConfManager : public QWidget
{
    Q_OBJECT
public:
    explicit ConfManager(kryomol::World* w,QWidget *parent = nullptr);


signals:
private:
    void InitTree();
private:
    kryomol::World* m_world;
    QTreeWidget* m_tree;
};

#endif // CONFMANAGER_H
