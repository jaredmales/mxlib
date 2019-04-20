# Definitions for building mxlib based applications
#
SELF_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
-include $(SELF_DIR)/../local/MxApp.mk
include $(SELF_DIR)/../mk/Common.mk

NEED_BLAS ?= yes
NEED_FFTW ?= yes
NEED_SOFA ?= yes
NEED_LEVMAR ?= yes
NEED_FITS ?= yes
NEED_BOOST ?= yes
NEED_GSL ?= yes
NEED_XPA ?= yes

# Provide default build options in case they weren't defined in local/MxApp.mk
ifeq ($UNAME,Darwin)  # macOS
    USE_BLAS_FROM ?= Accelerate
else
    USE_BLAS_FROM ?= mkl
endif

#default FFT is fftw
USE_FFT_FROM ?= fftw

# Configure includes and libraries based on build options
ifeq ($(USE_BLAS_FROM),mkl)
    BLAS_INCLUDES ?= -DMXLIB_MKL -m64 -I${MKLROOT}/include
    BLAS_LDFLAGS ?= -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed
    BLAS_LDLIBS ?= -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
endif
ifeq ($(USE_BLAS_FROM),ATLAS)
    #These are probably what you want for self-compiled atlas
    BLAS_INCLUDES ?= -I/usr/local/atlas/include
    BLAS_LDFLAGS ?= -L/usr/local/atlas/lib
    BLAS_LDLIBS ?= -llapack -lf77blas -lcblas -latlas -lgfortran
    
    #2018-08-13: These work for apt installed atlas on Ubuntu 18.04
    #BLAS_INCLUDES ?= -I/usr/include/x86_64-linux-gnu/
    #BLAS_LDFLAGS ?= -L/usr/lib/x86_64-linux-gnu/
    #BLAS_LDLIBS ?= -llapack -lf77blas -lcblas -latlas

endif
ifeq ($(USE_BLAS_FROM),Accelerate)
    BLAS_LDFLAGS ?= -framework Accelerate
endif

ifeq ($(USE_FFT_FROM),fftw)
    #Order matters, _threads first.
    FFT_LDLIBS += -lfftw3_threads -lfftw3f_threads -lfftw3l_threads -lfftw3 -lfftw3f  -lfftw3l 
endif

ifeq ($(NEED_SOFA),yes)
   SOFA_LIB = -lsofa_c
else
   SOFA_LIB =
endif

ifeq ($(NEED_LEVMAR),yes)
   LEVMAR_LIB = -llevmar 
else
   LEVMAR_LIB =
endif

ifeq ($(NEED_FITS),yes)
   FITS_LIB = -lcfitsio 
else
   FITS_LIB =
endif

ifeq ($(NEED_BOOST),yes)
   BOOST_LIB = -lboost_system -lboost_filesystem
else
   BOOST_LIB =
endif

ifeq ($(NEED_XPA),yes)
   XPA_LIB = -lxpa
else
   XPA_LIB = 
endif

ifeq ($(NEED_GSL),yes)
   GSL_LIB = -lgsl 
else
   GSL_LIB = 
endif

OPTIMIZE ?= -O3 -fopenmp -ffast-math
EXTRA_LDLIBS ?= $(SOFA_LIB) $(LEVMAR_LIB) $(FITS_LIB) $(BOOST_LIB) $(GSL_LIB) $(XPA_LIB)

ifneq ($(UNAME),Darwin)
    EXTRA_LDLIBS += -lrt
endif

EXTRA_LDFLAGS ?= -L$(PREFIX)/lib

ifeq ($(NEED_BLAS),yes)
    INCLUDES += $(BLAS_INCLUDES)
    EXTRA_LDLIBS += $(BLAS_LDLIBS)
    EXTRA_LDFLAGS += $(BLAS_LDFLAGS)
endif

ifeq ($(NEED_FFTW),yes)
   EXTRA_LDLIBS += $(FFT_LDLIBS)
   EXTRA_LDFLAGS += $(FFT_LDFLAGS)
endif

LDLIBS += $(EXTRA_LDLIBS) 
LDFLAGS += $(EXTRA_LDFLAGS)

CFLAGS += -std=c99 -fPIC $(INCLUDES) $(OPTIMIZE)
CXXFLAGS += -std=c++14 -fPIC $(INCLUDES) $(OPTIMIZE)

#This is needed to force use of g++ for linking
LINK.o = $(LINK.cc)

# Single-file app name can be supplied as `TARGET=`,
# or `t=` for short
TARGET ?= $(t)

all: $(TARGET) $(OTHER_OBJS)

$(TARGET):  $(TARGET).o  $(OTHER_OBJS)
	$(LINK.o)  -o $(TARGET) $(TARGET).o $(OTHER_OBJS) $(LDFLAGS) $(LDLIBS)
	
install: all
	install -d $(BIN_PATH)
	install $(TARGET) $(BIN_PATH)

.PHONY: clean
clean:
	rm -f $(TARGET)
	rm -f *.o 
	rm -f *~
