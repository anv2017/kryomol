#ifndef QJOBDYNWIDGET_H
#define QJOBDYNWIDGET_H

#include "qjobwidget.h"

class QConvWidget;
class QJobDynWidget : public QJobWidget
{
    Q_OBJECT
public:
    explicit QJobDynWidget(QWidget *parent = 0);
    void InitWidgets();


private:
    QConvWidget* m_dynwidget;
};

#endif // QJOBDYNWIDGET_H
