# Definitions for building mxlib based applications
#
SELF_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
-include $(SELF_DIR)/../local/MxApp.mk
include $(SELF_DIR)/../mk/Common.mk

# Provide default build options in case they weren't defined in local/MxApp.mk
ifeq ($UNAME,Darwin)  # macOS
    USE_BLAS_FROM ?= vecLib
else
    USE_BLAS_FROM ?= mkl
endif
USE_FFT_FROM ?= mkl

# Configure includes and libraries based on build options
ifeq ($(USE_BLAS_FROM),mkl)
    BLAS_INCLUDES ?= -m64 -I${MKLROOT}/include
    BLAS_LDFLAGS ?= -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed
    BLAS_LDLIBS ?= -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
endif
ifeq ($(USE_BLAS_FROM),ATLAS)
    BLAS_INCLUDES ?= -I/usr/local/atlas/include
    BLAS_LDFLAGS ?= -L/usr/local/atlas/lib
    BLAS_LDLIBS ?= -llapack -lf77blas -lcblas -latlas -lgfortran
endif

ifeq ($(USE_FFT_FROM),fftw)
    FFT_LDLIBS += -lfftw3f -lfftw3  -lfftw3l -lfftw3q -lfftw3f_threads -lfftw3_threads -lfftw3l_threads -lfftw3q_threads 
endif

OPTIMIZE ?= -O3 -fopenmp -ffast-math
EXTRA_LDLIBS ?= -lsofa_c -llevmar -lcfitsio -lrt -lboost_system -lboost_filesystem -lgsl
EXTRA_LDFLAGS ?= -L$(PREFIX)/lib

LDLIBS += $(BLAS_LDLIBS) $(FFT_LDLIBS) $(EXTRA_LDLIBS)
LDFLAGS += $(BLAS_LDFLAGS) $(FFT_LDFLAGS) $(EXTRA_LDFLAGS)

CFLAGS += -std=c99 -fPIC $(INCLUDES) $(OPTIMIZE)
CXXFLAGS += -std=c++14 -fPIC $(INCLUDES) $(OPTIMIZE)

# Single-file app name can be supplied as `TARGET=`,
# or `t=` for short
TARGET ?= $(t)

all: $(TARGET)

install: all
	install -d $(BIN_PATH)
	install $(TARGET) $(BIN_PATH)

.PHONY: clean
clean:
	rm $(TARGET)
	rm -f *.o 
	rm -f *~
