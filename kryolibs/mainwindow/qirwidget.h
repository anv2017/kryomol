/*****************************************************************************************
                            qirwidget.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QIRWIDGET_H
#define QIRWIDGET_H

#include <QWidget>
#include "ui_qirwidgetbase.h"
#include "kryomolcore_export.h"

class QJCDrawing;
class QPlotSpectrum;
class ConfManager;

namespace kryomol {
class World;
}

class KRYOMOLCORE_EXPORT QIRWidget : public QWidget, private Ui::QIRWidgetBase
{
    Q_OBJECT
public:
    QIRWidget(kryomol::World* w, QWidget* parent=0);

    ~QIRWidget();
    QPlotSpectrum* GetSpectrum() const { return m_spectrum; }
signals:
    void populations(const std::vector<double>& pops);
private:
    void InitButtons();
private slots:
    void OnShowConformers(bool );
    void OnGetPopulations();
private:
    QPlotSpectrum* m_spectrum;
    ConfManager* m_confmanager;
    kryomol::World* m_world;

};

#endif
