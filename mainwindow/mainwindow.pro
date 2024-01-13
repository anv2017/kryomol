#-------------------------------------------------
#
# Project created by QtCreator 2010-10-14T17:10:50
#
#-------------------------------------------------

QT += core gui opengl svg

TARGET = KryoMol
TEMPLATE = app

configfile = ../kryolibs/config/qryomol.pri
include($${configfile} )
CONFIG += opengl desktop

desktop{
    macx {
        DEFINES += D_MAC_WIN
    }
    win32 {
        DEFINES += D_MAC_WIN
    }
}

#unix {
#  !macx {
#  installdirfile=../installdir_unix.pri
#  include ($${installdirfile})
#  isEmpty(installdir) {
#  installdir=/usr/local/bin
#  }
#  target.path=$$installdir
#  INSTALLS+=target
#  }
#}


INCLUDEPATH += ../kryolibs/info ../kryolibs/core ../kryolibs/render ../kryolibs/plugin ../kryolibs/plugin/world ../kryolibs/plugin/tools ../kryolibs/cudatools ../kryolibs/tools ../kryolibs/parsers ../kryolibs/mainwindow  ../kryolibs/3dparty/qwt6/src ../kryolibs/3dparty/qwt6/designer ../kryolibs/3dparty/qwt6/lib ../kryolibs/wavelets

#detect shadow building and include folder with ui_*.h files and svninfo.h
!equals($$IN_PWD,$$OUT_PWD) {
message("shadow building detected")
INCLUDEPATH += $$OUT_PWD/../kryolibs/mainwindow  $$OUT_PWD/../kryolibs/info
}

win32{

  CONFIG(debug, debug|release){
    LIBS += -L../kryolibs/mainwindow/debug  -L../kryolibs/render/debug -L../kryolibs/plugin/debug \
            -L../kryolibs/parsers/debug -L../kryolibs/core/debug -L../kryolibs/cudatools/debug/ -L../kryolibs/tools/debug/ -L../kryolibs/wavelets/debug -L../kryolibs/3dparty/gl2ps/debug
  }
  CONFIG(release, debug|release){
    LIBS += -L../kryolibs/mainwindow/release  -L../kryolibs/render/release -L../kryolibs/plugin/release \
            -L../kryolibs/parsers/release -L../kryolibs/core/release -L../kryolibs/cudatools/release -L../kryolibs/tools/release -L../kryolibs/wavelets/release -L../kryolibs/3dparty/gl2ps/release
  }
}

LIBS += -L../kryolibs/mainwindow -L../kryolibs/script -L../kryolibs/render -L../kryolibs/plugin -L../kryolibs/parsers -L../kryolibs/core -L../kryolibs/cudatools -L../kryolibs/tools -L../kryolibs/wavelets -L../kryolibs/3dparty/gl2ps -L../kryolibs/3dparty/qwt/lib
LIBS += -lqryomolmainwindow -lqryomolrender -lqryomolplugin -lqryomolparsers -lqryomolcore  -lqryomoltools -lqryomolwavelets  -lqryomolgl2ps  #-llapack -lblas


#detect shadow building and include qwt6 folder
!equals($$IN_PWD,$$OUT_PWD) {
message("shadow building detected")
LIBS += $$OUT_PWD/../kryolibs/3dparty/qwt6/lib/libqwt.a
} else {
LIBS += ./kryolibs/3dparty/qwt6/lib/libqwt.a
}
win32 {
    LIBS += -lglu32 -lopengl32
}

unix {
!macx {
LIBS += -lGLU
}
}

cuda {
LIBS += -L../kryolibs/cudatools
LIBS += -lcudatools -lcudart
}


SOURCES += main.cpp \
    orcadialog.cpp \
    orcaengine.cpp \
    kryomolmainwindow.cpp \
    qjoboptwidget.cpp \
    qjobfreqwidget.cpp \
    qjobwidget.cpp \
    qjobuvwidget.cpp \
    qjobspwidget.cpp \
    qnicsgriddialog.cpp \
    qjobdynwidget.cpp

HEADERS  += kryomolmainwindow.h \
    orcadialog.h \
    orcaengine.h \
    qjoboptwidget.h \
    qjobfreqwidget.h \
    qjobwidget.h \
    qjobuvwidget.h \
    qjobspwidget.h \
    qnicsgriddialog.h \
    qjobdynwidget.h

RESOURCES += resources/qryomol.qrc

win32 {
RC_ICONS = resources/kryomol.ico
}

icons.files=resources/*.png
icons.path=mobkryolibs/icons

macx {
ICON = resources/kryomol.icns
}

libbasepath=
!equals($$IN_PWD,$$OUT_PWD) {
libbasepath=$$OUT_PWD/../kryolibs
} else {
libbasepath=$$IN_PWD/../kryolibs
}

unix {
POST_TARGETDEPS += $$libbasepath/tools/libqryomoltools.a \
                   $$libbasepath/core/libqryomolcore.a \
                   $$libbasepath/parsers/libqryomolparsers.a \
                   $$libbasepath/plugin/libqryomolplugin.a \
                   $$libbasepath/render/libqryomolrender.a \
                   $$libbasepath/mainwindow/libqryomolmainwindow.a \
}

win32 {
middlepath=
CONFIG(debug, debug|release) {
middlepath=Debug
}
CONFIG(release, debug|release){
middlepath=Release
}

POST_TARGETDEPS += $$libbasepath/tools/$$middlepath/libqryomoltools.a \
                   $$libbasepath/core/$$middlepath/libqryomolcore.a \
                   $$libbasepath/parsers/$$middlepath/libqryomolparsers.a \
                   $$libbasepath/plugin/$$middlepath/libqryomolplugin.a \
                   $$libbasepath/render/$$middlepath/libqryomolrender.a \
                   $$libbasepath/mainwindow/$$middlepath/libqryomolmainwindow.a \
}

FORMS += \
    qnicsgriddialog.ui



