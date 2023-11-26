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

class KRYOMOLCORE_EXPORT QIRWidget : public QWidget, private Ui::QIRWidgetBase
{
 Q_OBJECT
public:
    QIRWidget(QWidget* parent=0);

    ~QIRWidget();
    QPlotSpectrum* GetSpectrum() const { return m_spectrum; }
private:
  void InitButtons();
private:
  QPlotSpectrum* m_spectrum;

};

#endif
