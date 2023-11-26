#include "openingthread.h"
#include "mobdevwindow.h"

#include <QThread>
#include <QProgressDialog>

//OpeningThread::OpeningThread (QProgressDialog* progressDialog, QWidget* parent) : QThread (parent), m_progressDialog (progressDialog)
OpeningThread::OpeningThread (QWidget* parent) : QThread (parent)
{

}

OpeningThread::~OpeningThread()
{}

void OpeningThread::run()
{
    m_progressDialog = new QProgressDialog();
    m_progressDialog->setRange(0, 0);
    m_progressDialog->setValue(0);

    m_progressDialog->show();

    m_progressDialog->exec();

    //this->exec();
}

void OpeningThread::close()
{
    m_progressDialog->close();
}
