#ifndef ORCAENGINE_H
#define ORCAENGINE_H

#include <QTemporaryDir>

class QProcess;
class OrcaEngine : public QObject
{
    Q_OBJECT
public:
    OrcaEngine();
    void SetInputFile(const QString& string);
    void Run();
    //QString OutputFileName() const { return m_outputfilename; }
    QString OutputString() const { return m_outputstring; }
private slots:
    void OnJobFinished(int exitcode);
signals:
    void finished(int );
private:
    void Init();

private:
    QTemporaryDir* m_tmpdir;
    QDir m_jobdir;
    QString m_inputstring;
    QString m_outputstring;
    QProcess* m_process;
    QString m_jobname;
    QString m_outputfilename;

};

#endif // ORCAENGINE_H
