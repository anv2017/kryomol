/*****************************************************************************************
                            qdoubleeditbox.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QDOUBLEEDITBOX_H
#define QDOUBLEEDITBOX_H


#include <QWidget>

class QLineEdit;
class QPushButton;
class QDoubleEditBox : public QWidget
{
  Q_OBJECT
public:
    QDoubleEditBox(int first, int last, QWidget* parent=0);
    QDoubleEditBox(QWidget* parent=0);
    ~QDoubleEditBox();
signals:
  void limits(int,int);
public slots:
  void SetLimits(int, int);
private slots:
  void OnReturnPressed();
  void SetupConnections();
  void OnResetButton();
private:
  void SetValues(int first, int last);
  void InitWidget();
private:
 int m_first;
 int m_last;
 QLineEdit* m_toLineEdit;
 QLineEdit* m_fromLineEdit;
 QPushButton* m_resetButton;

};

#endif
