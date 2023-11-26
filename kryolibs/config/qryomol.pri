DEFINES += KRYOMOLSTATIC

QT += widgets
#CONFIG += cuda

CONFIG += timers
#CONFIG += noverbose
CONFIG += sse

macx {
CONFIG += veclib
}

cuda {
DEFINES += WITH_CUDA
CUDA_DIR = $$system(which nvcc | sed 's,/bin/nvcc$,,')
INCLUDEPATH += $$CUDA_DIR/include
LIBS += -L$$CUDA_DIR/lib64
}

win32 {
  mkl {
  INCLUDEPATH+=$$(MKLROOT)/include
  }
  CONFIG+= c++17
}


timers {
DEFINES += WITH_TIMERS
}

noverbose {
DEFINES += QT_NO_DEBUG_OUTPUT
}

sse {
 unix {
 QMAKE_CXXFLAGS += -mmmx -mfpmath=sse -msse -msse2 -msse3  -Wno-narrowing
 }
 win32 {
 QMAKE_CXXFLAGS += -mmmx -mfpmath=sse -msse -msse2 -msse3  -Wno-narrowing
 }
 DEFINES += WITH_SIMD
}

macx {
# altivec {
# QMAKE_CXXFLAGS += -faltivec
# }
}

