TEMPLATE = subdirs

#include first the cudatools diretory otherwise it will complain.
SUBDIRS +=  3dparty tools #info

cuda{
    SUBDIRS += cudatools
}

SUBDIRS += core render parsers plugin mainwindow script wavelets
