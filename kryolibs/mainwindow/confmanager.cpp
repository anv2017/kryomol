#include "confmanager.h"

#include <QTreeWidget>
#include <QVBoxLayout>
#include <QFrame>
#include <QLabel>
#include <QCheckBox>

#include "world.h"
#include "molecule.h"

ConfManager::ConfManager(kryomol::World* w,QWidget* parent)
    : QWidget{parent}, m_world(w)
{
    QVBoxLayout* vbox = new QVBoxLayout(this);
    m_tree= new QTreeWidget(this);
    vbox->addWidget(m_tree);
    InitTree();
}

void ConfManager::Refresh()
{
    InitTree();
}

void ConfManager::InitTree()
{
    m_tree->clear();
    QStringList headers;
    headers << "" << "Color" << "Population" << "Include" << "Show";
    m_tree->setColumnCount(headers.size());
    m_tree->setHeaderLabels(headers);

    //Do not make this item a child of QTreeWidget or it will be used as
    //a normal QTreeWidget Item
    QTreeWidgetItem* spitem = new QTreeWidgetItem(m_tree);
    spitem->setText(0,"Spectra");
    spitem->setExpanded(true);

    /*QCheckBox* sp3 = new QCheckBox(m_tree);
    QCheckBox* sp4 = new QCheckBox(m_tree);
    sp3->setChecked(true);
    sp4->setChecked(true);
    m_tree->setItemWidget(spitem,3, sp3);
    m_tree->setItemWidget(spitem,4, sp4);*/
    spitem->setFlags(spitem->flags() | Qt::ItemIsAutoTristate);
    spitem->setCheckState(3,Qt::Checked);


    std::vector<QTreeWidgetItem*> items;
    QTreeWidgetItem* avitem = new QTreeWidgetItem(spitem);
    items.push_back(avitem);
    avitem->setText(0,"Average");
    QFrame* avfr = new QFrame(m_tree);
    avfr->setAutoFillBackground(true);
    avfr->setStyleSheet( QString("background-color: black"));
    spitem->addChild(avitem);
    m_tree->setItemWidget(avitem,1,avfr);


    for(size_t fidx=0;fidx<m_world->CurrentMolecule()->Frames().size();++fidx)
    {
        const kryomol::Frame& frame=m_world->CurrentMolecule()->Frames()[fidx];
        const double& pop=m_world->CurrentMolecule()->Populations()[fidx];
        QTreeWidgetItem* fitem = new QTreeWidgetItem(spitem);
        spitem->addChild(fitem);
        items.push_back(fitem);
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
        fitem->setText(2,QString::number(pop,'g',3));

    }

    for(auto p : items)
    {
        p->setCheckState(3,Qt::Checked);
        p->setCheckState(4,Qt::Checked);
    }

    for(int i=0;i<m_tree->columnCount();++i)
    {
        m_tree->resizeColumnToContents(i);
    }

    connect(m_tree,&QTreeWidget::itemChanged,this,&ConfManager::OnItemChanged);
}

void ConfManager::OnItemChanged(QTreeWidgetItem* item)
{
    qDebug() << "item here" << item << endl;
    if ( item == m_tree->topLevelItem(0) ) return;
    bool b;
    std::vector<bool> vb;

    for(int i=0;i<m_tree->topLevelItem(0)->childCount();++i)
    {
        QTreeWidgetItem* p=m_tree->topLevelItem(0)->child(i);
        if ( i == 0 )
        {
            b= ( p->checkState(4) == Qt::Checked );
        }
        else
        {
            vb.push_back( p->checkState(4) == Qt::Checked );
        }
    }

    std::cout << "b" << b << std::endl;
    for(auto x : vb)
    {
        std::cout << x << std::endl;
    }

    emit visible(b,vb);
}

