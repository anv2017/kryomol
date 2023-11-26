/** code inspired on qtcreator.blogspot.com */


#include <iostream>
#include <QProcess>
#include <QFile>
#include <QTextStream>

int main(int argc,char** argv)
{
    if ( argc != 2)
    {
        std::cout << "Usage: getsvninfo directory" << std::endl;
        return 1;
    }

    QFile file("svnversion.h");
    if (!file.open(QIODevice::WriteOnly))

    {
        std::cout << "Unable to open svnversion.h file" << std::endl;
        return 1;
    }

    QProcess process;
    std::cout << argv[1] << std::endl;
    process.start("svnversion", QStringList() << "-n" << argv[1]);

    if ( process.waitForStarted() && process.waitForReadyRead())
    {
        QTextStream voutput(&file);
        voutput << "#define " << "QRYOMOLSVNVERSION " << '"' << process.readAllStandardOutput().constData() << '"' << endl;
    }

    return 0;
}
