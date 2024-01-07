/*****************************************************************************************
                            main.cpp  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/


#include "kryomolmainwindow.h"
#include "qryomolapp.h"
#include "iostream"


int main(int argc, char *argv[])
{
    kryomol::KryoMolApplication a(argc, argv);

    KryoMolMainWindow *mw=NULL;
    int nargs;

    nargs=argc-1;

    if(nargs > 1)
    {
        std::cerr << "Excessive Number of arguments" << std::endl;
        exit(1);
    }

    if(nargs == 1)  //several files
    {

        mw= new KryoMolMainWindow();
        a.SetFile(argv[1]);

    }
    else
        mw = new KryoMolMainWindow();

    mw->showMaximized();

    return a.exec();

}

