#include "confmanager.h"

#include <QTreeWidget>
#include <QVBoxLayout>
#include <QFrame>
#include <QLabel>

#include "world.h"
#include "molecule.h"

ConfManager::ConfManager(kryomol::World* w,QWidget* parent)
    : m_world(w), QWidget{parent}
{
    QVBoxLayout* vbox = new QVBoxLayout(this);
    m_tree= new QTreeWidget(this);
    vbox->addWidget(m_tree);
    InitTree();
}

void ConfManager::InitTree()
{
    QStringList headers;
    headers << "Mol" << "Color" << "Population" << "Include" << "Show";
    m_tree->setColumnCount(headers.size());
    m_tree->setHeaderLabels(headers);
    QTreeWidgetItem* witem = new QTreeWidgetItem(m_tree);
    witem->setText(0,"Average");

    m_tree->addTopLevelItem(witem);

    for(size_t fidx=0;fidx<m_world->CurrentMolecule()->Frames().size();++fidx)
    {
        const kryomol::Frame& frame=m_world->CurrentMolecule()->Frames()[fidx];
        QTreeWidgetItem* fitem = new QTreeWidgetItem(witem);
        fitem->setText(0,QString::number(fidx+1));
        QFrame* fr = new QFrame(m_tree);
        fr->setFrameShape(QFrame::Box);
        fr->setAutoFillBackground(true);

        float hf,sf,lf;
        frame.GetColor(hf,sf,lf);
        int h=qRound(hf*360);
        int s=qRound(sf*100);
        int l=qRound(lf*100);

        qDebug() <<  QString("background-color: hsl(%1, %2%, %3%)").arg(hf).arg(sf).arg(lf) << endl;

        fr->setStyleSheet( QString("background-color: hsl(%1, %2%, %3%)").arg(h).arg(s).arg(l));
        m_tree->setItemWidget(fitem,1,fr);
        witem->addChild(fitem);



    }

    for(int i=0;i<m_tree->columnCount();++i)
    {
        m_tree->resizeColumnToContents(i);
    }
}


