TEMPLATE=app
CONFIG += qt

TARGET=getsvninfo

SOURCES += getsvninfo.cpp

HEADERS += qryomolinfo.h

#build console application on macosx
macx:CONFIG -=app_bundle


#detect shadow building and include output folder
!equals($$IN_PWD,$$OUT_PWD) {
message("shadow building detected in info")
INCLUDEPATH += $$OUT_PWD
message($$IN_PWD)
message($$OUT_PWD)


#create an empty svnversion.h
QMAKE_EXTRA_TARGETS=version
PRE_TARGETDEPS+=version.h
version.target=version.h
win32 {
version.commands = echo 2>svnversion.h
}
unix {
version.commands = touch svnversion.h
}

unix:QMAKE_POST_LINK=./$$TARGET $$IN_PWD
win32 {
CONFIG(debug, debug|release){
QMAKE_POST_LINK=$$OUT_PWD/Debug/getsvninfo.exe $$IN_PWD
}
CONFIG(release, debug|release){
QMAKE_POST_LINK=$$OUT_PWD/Release/getsvninfo.exe $$IN_PWD

}
}

} else {
unix:QMAKE_POST_LINK=./$$TARGET .

win32 {
        CONFIG(debug, debug|release){
        QMAKE_POST_LINK=Debug/$$(TARGET).exe ..
        }
        CONFIG(release, debug|release){
        QMAKE_POST_LINK=Release/$$(TARGET).exe ..
        }
}

}
