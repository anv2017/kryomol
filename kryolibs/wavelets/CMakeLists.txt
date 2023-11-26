
include_directories(${CMAKE_SOURCE_DIR}/src/kryonmr/kryonmr/kryonmrparsers ${CMAKE_SOURCE_DIR}/src/tools ${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(singleprecisionwlt_STAT_SRCS
    daubechies.c
    dwt.c
    error.c
    stream.c
    wavelet.c)

kde4_add_library(singleprecisionwlt STATIC ${singleprecisionwlt_STAT_SRCS})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")

########### install files ###############




#original Makefile.am contents follow:

#
#noinst_LIBRARIES = libsingleprecisionwlt.a
#
#INCLUDES = -I$(top_srcdir)/src/kryonmr/kryonmr/kryonmrparsers 	-I$(top_srcdir)/src/tools $(all_includes)
#
#libsingleprecisionwlt_a_CFLAGS = -fPIC
#
#noinst_HEADERS = gsl_errno.h gsl_message.h gsl_types.h gsl_wavelet.h
#
#libsingleprecisionwlt_a_SOURCES = daubechies.c dwt.c error.c stream.c wavelet.c
