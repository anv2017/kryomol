#ifndef VERSION_H
#define VERSION_H

#include <QString>
#include <QStringList>
#include <QTextStream>

#ifndef KRYOMOLSVNVERSION
#define KRYOMOLSVNVERSION ""
#endif

#define KRYOMOLMAJORVERSION 2
#define KRYOMOLMINORVERSION 0
#define KRYOMOLSUBMINORVERSION "2024.1"

namespace kryomol
{
class KryoMolInfo
{
public:
    static QString Version()
    {
        QString str;
        QTextStream stream(&str);
        stream << KRYOMOLMAJORVERSION << "." << KRYOMOLMINORVERSION << "-" << KRYOMOLSUBMINORVERSION;
        return str;
    }

    static QString SvnVersion()
    {
        return KRYOMOLSVNVERSION;
    }

    static QStringList Authors()
    {
        QStringList authorlist;
        authorlist << "Armando Navarro-Vázquez";
        return authorlist;

    }

    static QStringList Collaborators()
    {
        QStringList authorlist;
        authorlist << "Noa Campos López";
        return authorlist;
    }

    static QString Company()
    {
        return "";
    }

    static QStringList ThirdPartyPackages()
    {
        QStringList tparty;
        tparty << QString("Qwt6: http://qwt.sourceforge.net/");
        tparty << QString("GL2PS: http://www.geuz.org/gl2ps/");
        tparty << QString("SSE math library: http://gruntthepeon.free.fr/ssemath/");
        tparty << QString("GNU gsl library: https://www.gnu.org/software/gsl/");
        return tparty;
    }

    static QString Copyright()
    {
        return "This program is free software. It can be redistribute under the terms of the GPL v2 or any later version";

    }

    static QString Disclaimer()
    {
        return "The program is provided AS IS with NO WARRANTY OF ANY KIND,"
               "INCLUDING THE WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.";
    }
};

}
#endif
