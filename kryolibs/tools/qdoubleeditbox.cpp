/*****************************************************************************************
                            qdoubleeditbox.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <QLineEdit>
#include <QPushButton>
#include <QLayout>
#include <QLabel>

#include "qdoubleeditbox.h"

QDoubleEditBox::QDoubleEditBox(int first, int last, QWidget* parent )
 : QWidget(parent),  m_first(first), m_last(last)
{

  InitWidget();
  SetValues(first,last);
  SetupConnections();
}

QDoubleEditBox::QDoubleEditBox(QWidget* parent )
 : QWidget(parent)
{
   InitWidget();
   SetupConnections();
}

void QDoubleEditBox::InitWidget()
{
  QHBoxLayout* box= new QHBoxLayout(this);
  QLabel* fromlabel= new QLabel("From",this);
  QLabel* tolabel = new QLabel("to",this);
  m_fromLineEdit= new QLineEdit(this);
  m_toLineEdit= new QLineEdit(this);
  m_resetButton = new QPushButton(this);
  m_resetButton->setText("Reset");
  box->addWidget(fromlabel);
  box->addWidget(m_fromLineEdit);
  box->addWidget(tolabel);
  box->addWidget(m_toLineEdit);
  box->addWidget(m_resetButton);
}

QDoubleEditBox::~QDoubleEditBox()
{
}

void QDoubleEditBox::SetLimits(int first,int last)
{
  m_first=first;
  m_last=last;
  SetValues(first,last);
}
void QDoubleEditBox::OnReturnPressed()
{
  int first=m_fromLineEdit->text().toInt();
  int last=m_toLineEdit->text().toInt();
  emit limits(first,last);
}

void QDoubleEditBox::OnResetButton()
{
  SetValues(m_first,m_last);
  emit limits(m_first,m_last);
}

void QDoubleEditBox::SetValues(int first, int last)
{
   QString snum;
   snum.setNum(first);
   m_fromLineEdit->setText(snum);
   snum.setNum(last);
   m_toLineEdit->setText(snum);
}

void QDoubleEditBox::SetupConnections()
{
  connect(m_fromLineEdit,SIGNAL(returnPressed()),this,SLOT(OnReturnPressed()));
  connect(m_toLineEdit,SIGNAL(returnPressed()),this,SLOT(OnReturnPressed()));
  connect(m_resetButton,SIGNAL(clicked()),this,SLOT(OnResetButton()));
}

