#ifndef OPENINGTHREAD_H
#define OPENINGTHREAD_H

#include <QWidget>
#include <QThread>

class QProgressDialog;

class OpeningThread : public QThread
{
    Q_OBJECT

public:
    OpeningThread(QWidget* parent = 0);//(QProgressDialog* progressDialog, QWidget* parent = 0);
    ~OpeningThread();

    void close();

protected:
    void run();

private:
    QProgressDialog* m_progressDialog;
};

#endif // OPENINGTHREAD_H
