#include "orcaengine.h"

#include <QTemporaryDir>
#include <QProcess>
#include <QTextStream>
#include <QDebug>

OrcaEngine::OrcaEngine()
{
    // m_jobdir=QDir("/Users/armando/orca_test");
    // m_jobname="myfirstjob";
}

void OrcaEngine::SetInputFile(const QString& string)
{
    m_inputstring=string;
}

void OrcaEngine::Run()
{
    m_tmpdir = new QTemporaryDir();

    if ( !m_tmpdir->isValid() ) throw QString("Invalid temporary directory");

    QString filepath=m_tmpdir->path()+QDir::separator()+"INPUT";

    qDebug() << "filepath=" << filepath;

    QFile file(filepath);

    if ( ! file.open(QFile::WriteOnly) ) throw QString("Cannot open file for writing");

    QTextStream stream(&file);

    stream << m_inputstring;

    m_process = new QProcess();
    m_process->setWorkingDirectory(m_tmpdir->path());
    connect(m_process,SIGNAL(finished(int)),this,SLOT(OnJobFinished(int )));
    QStringList arguments;
    arguments << "INPUT";
    QString orcapath;
#ifdef Q_OS_WIN
    orcapath="c:/Orca/orca";
#else
#ifdef Q_OS_MAC
    orcapath="/Users/armando/Library/Orca421/orca";
#endif
#endif
    m_process->start(orcapath,arguments);

    if ( m_process->waitForStarted() )
    {
        qDebug() << "job started succesfully";
    }

}

void OrcaEngine::OnJobFinished(int exitcode)
{
    m_outputstring=m_process->readAllStandardOutput();
    m_outputfilename=m_jobdir.absolutePath()+QDir::separator()+m_jobname+".out";
    QFile outputfile(m_outputfilename);
    outputfile.open(QFile::WriteOnly);
    QTextStream outstream(&outputfile);
    outstream << m_outputstring;
    outputfile.close();
    emit finished(exitcode);
}
