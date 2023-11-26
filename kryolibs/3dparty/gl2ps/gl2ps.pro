configfile = ../../config/qryomol.pri
include($${configfile} )
TEMPLATE = lib
CONFIG +=  static opengl
TARGET = qryomolgl2ps
HEADERS += gl2ps.h
SOURCES += gl2ps.c
#DEFINES += GL2PSDLL GL2PSDLL_EXPORTS
!macx {
target.path = "$$(KRYOMOL_DIR)/lib"
} else {
QMAKE_LFLAGS_SONAME  = -Wl,-install_name,@executable_path/../lib/
target.path = "$$(KRYOMOL_DIR)/$${APPNAME}.app/Contents/lib"
}
#INSTALLS += target
win32 {
DLLDESTDIR = "$$(KRYOMOL_DIR)/bin"
}
