/*****************************************************************************************
                            qdoubleslider.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QDOUBLESLIDER_H
#define QDOUBLESLIDER_H

#include <QSlider>
#include <QWidget>

class QDoubleSlider : public QSlider
{
    Q_OBJECT
public:
    QDoubleSlider(QWidget* parent){}
    ~QDoubleSlider(){}

    void setDoubleMaximum(double max) {m_max = max;}
    void setDoubleMinimum(double min) {m_min = min;}

public slots:
    void OnValueChanged(int value);
    void OnSliderReleased();

signals:
    void valueChanged(double v);    

private:
    double m_min;
    double m_max;
    double m_value;
    double m_singlestep;
    int m_divisions;
};

#endif // QDOUBLESLIDER_H
