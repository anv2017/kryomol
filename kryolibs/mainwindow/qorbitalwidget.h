/*****************************************************************************************
                            qorbitalwidget.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QORBITALWIDGET_H
#define QORBITALWIDGET_H

#include <QWidget>

#include "ui_qorbitalwidgetbase.h"
#include "renderorbitals.h"
#include "density.h"
#include "qplotspectrogram.h"
#include "world.h"


class QButtonGroup;

class QOrbitalWidget : public QWidget, private Ui::QOrbitalWidgetBase
{
    Q_OBJECT

public:
    explicit QOrbitalWidget(QWidget *parent = 0);
    ~QOrbitalWidget();
    void SetWorld(kryomol::World* w);
    void Init();
    void ListOrbitals();
    void SetHomo(int h) {m_homo=h;}
    void SetLumo(int l) {m_lumo=l;}
    void SetBeta(bool b) {m_beta=b;}
    int Orbital() {return m_orbital;}
    void SetRenderOrbitals(kryomol::RenderOrbitals render);
    void HideAllButtons();
    kryomol::Density GetDensity() {return m_render.DensityData();}

signals:
    void drawDensity(bool );
    void offtransitions(bool );
    void transparenceChange(float );
public slots:
    void OnSetFrame(size_t frame);
private slots:
    void OnIsovalueSliderChange(double );
    void OnTransparenceSliderChange(int );
    void OnOrbitalChange(QTreeWidgetItem* );
    void OnBetaOrbitalChange(QTreeWidgetItem* );
    void OnThresholdAOSelectorChange(double );
    void OnGridResolutionChange(int );
    void OnShowTransitionChange(int );
    void OnShowDensityChange(int );
    void OnShowHomo();
    void OnShowLumo();
    void OnShowTotalDensity();
    void OnShowSpinDensity();
    void OnShowContoursFrame();
    void OnChangeAxisPlane(bool );
    void OnChangeAxisValue(double );
    void OffButtons(bool );

private:    
    bool m_beta;
    bool m_bshowhomo;
    bool m_bshowlumo;
    bool m_bshowtotaldensity;
    bool m_bshowspindensity;
    bool m_bonshow;
    bool m_bshowcontours;
    int m_homo;
    int m_lumo;
    int m_orbital;
    int m_betaorbital;
    float m_gridresolution;
    float m_isovalue;
    double m_threshold;
    kryomol::RenderOrbitals m_render;
    QPlotSpectrogram *m_plotspectrogram;
    kryomol::World* m_world;

};

#endif // QORBITALWIDGET_H
