/*****************************************************************************************
                            qmolecularlistcontrol.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QMOLECULARLISTCONTROL_H
#define QMOLECULARLISTCONTROL_H

#include "ui_qmolecularlistcontrolbase.h"

namespace kryomol
{
class World;
}
class QMolecularTreeWidget;
class QButtonGroup;

class QMolecularListControl : public QWidget , private Ui::QMolecularListControlBase
{
Q_OBJECT
public:
    QMolecularListControl(QWidget* parent =0);

    ~QMolecularListControl();
    void SetWorld(kryomol::World* world);
    void Init();
private slots:
    void OnChangeEnsamble(int id);
    void OnChangeTemperature();
    void OnCalculateButton();
    void OnReadPopulationsButton();
    void OnReadEnergiesButton();
    void OnPastePopulationsButton();
    void OnPasteEnergiesButton();
private:
    void ReadEnergiesFromStream(std::istream& stream);
    void ReadPopulationsFromStream(std::istream& stream);
private:
    QMolecularTreeWidget* m_listview;
    kryomol::World* m_world;
    QButtonGroup* m_ensamblegroup;
};

#endif
