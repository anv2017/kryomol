#ifndef QBASEJOBWIDGET_H
#define QBASEJOBWIDGET_H


#include <QWidget>
#include <QMainWindow>


namespace nmrdev
{
class World;
class NMRDevApplication;
class GLVisor;
};

class QBaseJobWidget : public QWidget
{
public:
    QBaseJobWidget(QWidget* parent=0);
    QBaseJobWidget(const char* file, QWidget* parent=0);
    ~QBaseJobWidget();

private:
    void Init();

protected:
    /** absolute path to the opened file*/
    QString m_file;
    nmrdev::NMRDevApplication* m_app;
    /** pointer to the simulation world*/
    nmrdev::World* m_world;
    nmrdev::GLVisor* m_visor;


};

#endif // QBASEJOBWIDGET_H
