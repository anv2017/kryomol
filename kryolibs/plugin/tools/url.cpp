/*****************************************************************************************
                            url.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <QUrl>
#include <QProcess>
#include <QMessageBox>
#ifdef Q_OS_MACX
#include <Carbon/Carbon.h>
#endif // Q_OS_MACX

#if QT_VERSION >= 0x040200
#include <QDesktopServices>
#endif

#include "url.h"
#include <iostream>

void kryomol::OpenURL(const QString& url)
{

#if QT_VERSION <= 0x040200
 QString myurl=url;
#ifdef Q_OS_MACX //from the Qt mailing list
/*
* Use Internet Config to hand the URL to the appropriate applicaton.
* as set by the user in the Internet Preferences pane.
*/
ICInstance icInstance;
uint psiSignature = 'Psi ';
OSStatus error = ::ICStart ( &icInstance, psiSignature );
if ( !error ) 
{
ConstStr255Param hint ( 0x0 );
const char *data = url.toLatin1(); 
long length = url.length();
long start = ( 0 );
long end ( length );
// Don't bother testing return value ( error ); Launched application
// will report problems.
::ICLaunchURL ( icInstance, hint, data, length, &start, &end );
ICStop ( icInstance );
}
return;
#endif // Q_OS_MACX

#ifdef Q_WS_X11
  QProcess* pr= new QProcess();
  char* browser=getenv("BROWSER");
  QString sbrowser(browser);
  QStringList browserlist;
  std::cout << "browser" << browser << std::endl;
  if ( browser == 0 )
  {
    //well other browsers

    browserlist << "firefox" << "mozilla-firefox" << "opera" << "mozilla" << "netscape" << "konqueror";
  }
  else
  {
    browserlist = sbrowser.split(":");
  }

  QStringList::Iterator it;
  QStringList args;
  args << myurl;
  for(it=browserlist.begin();it!=browserlist.end();it++)
  {

    pr->start(*it,args);// return;
    if ( pr->waitForStarted() )
      return;
  }
#endif

#ifdef Q_WS_WIN
  
 ShellExecute( 0,L"open",(TCHAR*)url.unicode(),NULL,NULL,SW_SHOW);
 return;

#endif
  QMessageBox mb;
  mb.setText("Sorry cannot find any browser on your system");
  mb.exec();
#else
    QDesktopServices::openUrl(QUrl(url));
#endif

}
