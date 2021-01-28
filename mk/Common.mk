# Common path and make variable definitions
#
# NOTE: This file should only be edited in mxlib/local, not in the root mxlib directory.
#
SELF_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
-include $(SELF_DIR)/../local/Common.mk


# Set these to no in local/Common.mk if you never need them
# or in your local Makefile
NEED_CUDA ?= yes


UNAME ?= $(shell uname)
ifeq ($(UNAME),Darwin)
	CFLAGS += -D_BSD_SOURCE
	CXXFLAGS += -D_BSD_SOURCE
else
	CFLAGS += -D_XOPEN_SOURCE=700
	CXXFLAGS += -D_XOPEN_SOURCE=700
endif
PREFIX ?= $(HOME)
BIN_PATH ?= $(PREFIX)/bin
LIB_PATH ?= $(PREFIX)/lib
INCLUDE_PATH ?= $(PREFIX)/include
LIB_SOFA ?= $(LIB_PATH)/libsofa_c.a
ARFLAGS ?= rvs

INCLUDES += -I$(INCLUDE_PATH)

OPTIMIZE ?= -O3 -fopenmp -ffast-math

CFLAGS += -std=c99 -fPIC
CXXFLAGS += -std=c++14 -fPIC

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
    #FFT_LDLIBS ?= -lfftw3_threads -lfftw3f_threads -lfftw3l_threads -lfftw3 -lfftw3f  -lfftw3l 
    
    #with new combined-threads:
    FFT_LDLIBS ?= -lfftw3 -lfftw3f -lfftw3l -lfftw3q 
endif

FITS_LIB = -lcfitsio 

BOOST_LIB = -lboost_system -lboost_filesystem

GSL_LIB = -lgsl 

SOFA_LIB = -lsofa_c

EXTRA_LDLIBS ?= $(FITS_LIB) $(BOOST_LIB) $(GSL_LIB) $(SOFA_LIB)

ifneq ($(UNAME),Darwin)
    EXTRA_LDLIBS += -lrt
endif

EXTRA_LDFLAGS ?= -L$(PREFIX)/lib

#BLAS:
INCLUDES += $(BLAS_INCLUDES)
EXTRA_LDLIBS += $(BLAS_LDLIBS)
EXTRA_LDFLAGS += $(BLAS_LDFLAGS)

#FFTW:
EXTRA_LDLIBS += $(FFT_LDLIBS)
EXTRA_LDFLAGS += $(FFT_LDFLAGS)


LDLIBS += $(EXTRA_LDLIBS) 
LDFLAGS += $(EXTRA_LDFLAGS)

CFLAGS += $(INCLUDES) $(OPTIMIZE)
CXXFLAGS += $(INCLUDES) $(OPTIMIZE)

#This is needed to force use of g++ for linking
LINK.o = $(LINK.cc)

ifeq ($(NEED_CUDA),yes)

   CXXFLAGS += -DEIGEN_NO_CUDA
   
   HOST_ARCH   := $(shell uname -m)
   CUDA_TARGET_ARCH = $(HOST_ARCH)
   ifneq (,$(filter $(CUDA_TARGET_ARCH),x86_64 aarch64 ppc64le armv7l))
       ifneq ($(CUDA_TARGET_ARCH),$(HOST_ARCH))
           ifneq (,$(filter $(CUDA_TARGET_ARCH),x86_64 aarch64 ppc64le))
               TARGET_SIZE := 64
           else ifneq (,$(filter $(CUDA_TARGET_ARCH),armv7l))
               TARGET_SIZE := 32
           endif
       else
           TARGET_SIZE := $(shell getconf LONG_BIT)
       endif
   else
       $(error ERROR - unsupported value $(CUDA_TARGET_ARCH) for TARGET_ARCH!)
   endif
   
   # operating system
   HOST_OS   := $(shell uname -s 2>/dev/null | tr "[:upper:]" "[:lower:]")
   TARGET_OS ?= $(HOST_OS)
   ifeq (,$(filter $(TARGET_OS),linux darwin qnx android))
       $(error ERROR - unsupported value $(TARGET_OS) for TARGET_OS!)
   endif
   
   # host compiler
   ifeq ($(TARGET_OS),darwin)
       ifeq ($(shell expr `xcodebuild -version | grep -i xcode | awk '{print $$2}' | cut -d'.' -f1` \>= 5),1)
           HOST_COMPILER ?= clang++
       endif
   endif

   HOST_COMPILER ?= g++
   NVCC          := nvcc -ccbin $(HOST_COMPILER)
   
   # internal flags
   NVCCFLAGS   := -m${TARGET_SIZE}
   NVCCFLAGS   +=  -DEIGEN_NO_CUDA -DMXLIB_MKL
   
   # build flags
   ifeq ($(TARGET_OS),darwin)
       LDFLAGS += -rpath $(CUDA_PATH)/lib
       CXXFLAGS += -arch $(HOST_ARCH)
   endif
   
   
   # Debug build flags
   ifeq ($(dbg),1)
         NVCCFLAGS += -g 
         BUILD_TYPE := debug
   else
         BUILD_TYPE := release
   endif
   
   ALL_CCFLAGS :=
   ALL_CCFLAGS += $(NVCCFLAGS)
   ALL_CCFLAGS += $(EXTRA_NVCCFLAGS)
   ALL_CCFLAGS += $(addprefix -Xcompiler ,$(CXXFLAGS))
   ALL_CCFLAGS += $(addprefix -Xcompiler ,$(EXTRA_CCFLAGS))
   ALL_CCFLAGS +=
   
   ALL_LDFLAGS :=
   ALL_LDFLAGS += $(ALL_CCFLAGS)
   ALL_LDFLAGS += $(addprefix -Xlinker ,$(LDFLAGS))
   ALL_LDFLAGS += $(addprefix -Xlinker ,$(LDLIBS))
   

   #build any cu and cpp files through NVCC as needed
   %.o : %.cu
	$(EXEC) $(NVCC) $(ALL_CCFLAGS) $< -c -o $@

#   %.o : %.cpp
#	$(EXEC) $(NVCC) $(ALL_CCFLAGS) $< -c -o $@


   #Finally we define the cuda libs for linking
   CUDA_LIBS ?= -L/usr/local/cuda/lib64/ -lcudart -lcublas -lcufft -lcurand
   
   
else
   CUDA_LIBS ?= 
endif
