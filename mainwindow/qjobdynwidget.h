#ifndef QJOBDYNWIDGET_H
#define QJOBDYNWIDGET_H

#include "qjobwidget.h"

class QConvWidget;
class QJobDynWidget : public QJobWidget
{
    Q_OBJECT
public:
    explicit QJobDynWidget(QWidget *parent = nullptr);
    void InitWidgets();


private:
    QConvWidget* m_dynwidget;
};

#endif // QJOBDYNWIDGET_H
