# Common path and make variable definitions
#
# NOTE: This file should only be edited in mxlib/local, not in the root mxlib directory.
#
SELF_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
-include $(SELF_DIR)/../local/Common.mk

UNAME ?= $(shell uname)
# environments setting -Wl,-dead_strip_dylibs by default break MKL link options
# so we remove it by replacing it with zero if it's there; later parts of the
# Makefile then add to this adjusted LDFLAGS value
LDFLAGS := $(LDFLAGS:-Wl,-dead_strip_dylibs=)

# Set to no in local/Common.mk or in your local Makefile if you never need CUDA
NEED_CUDA ?= yes

# Set to yes in local/Common.mk or in your local Makefile if you need MILK functions
NEED_MILK ?= no

ifeq ($(UNAME),Darwin)
	CFLAGS += -D_BSD_SOURCE
	CXXFLAGS += -D_BSD_SOURCE
    # CUDA's effectively unsupported on macOS, so don't bother:
    NEED_CUDA = no
    LIB_SUFFIX = dylib
else
	CFLAGS += -D_XOPEN_SOURCE=700
	CXXFLAGS += -D_XOPEN_SOURCE=700
    LIB_SUFFIX = so
    LDFLAGS += -Wl,-rpath=$(PREFIX)/lib
endif
PREFIX ?= $(HOME)
BIN_PATH ?= $(PREFIX)/bin
LIB_PATH ?= $(PREFIX)/lib
INCLUDE_PATH ?= $(PREFIX)/include
INCLUDES += -I$(INCLUDE_PATH)
ARFLAGS ?= rvs

DEFAULT_OPTIMIZATIONS = -O3 -ffast-math
ifeq ($(UNAME),Linux)
    DEFAULT_OPTIMIZATIONS += -fopenmp
endif
OPTIMIZE ?= $(DEFAULT_OPTIMIZATIONS)

CFLAGS += -std=c99 -fPIC

#Set to -std=c++14 if 17 not supported
CXXVERSION ?= -std=c++17

CXXFLAGS += $(CXXVERSION) -fPIC

USE_FFT_FROM ?= fftw
USE_BLAS_FROM ?= openblas

#***************************************************************#
# Configure includes and libraries based on build options
#***************************************************************#


################################################
## BLAS and Lapack
###############################################

# MKL
ifeq ($(USE_BLAS_FROM),mkl)
    $(if ${MKLROOT},,$(warning No value set for environment variable $$MKLROOT))
    BLAS_INCLUDES ?= -DMXLIB_MKL -m64 -I${MKLROOT}/include
    ifeq ($(UNAME),Darwin)
        BLAS_LDFLAGS ?= -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib
        BLAS_LDLIBS ?= -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
    endif
    ifeq ($(UNAME),Linux)
        BLAS_LDFLAGS ?= -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed
        BLAS_LDLIBS ?= -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
    endif
endif

# Open BLAS
# Provides several ways to compile against openblas.  For a system-wide blas install (i.e. the one
# chosen with update-alternatives) USE_BLAS_FROM can just be `openblas` or `openblas64`.
#
# A specific openblas version can be chosen on Ubuntu with the choices openblas64-serial/openmp/pthread
# or openblas-serial/openmp/pthread.  The defaults here are specified for Ubuntu 22.04.  Will probably
# work on any Ubuntu.  This can be fully customized with OPENBLAS_PKG_CONFIG_NAME, OPENBLAS_PKG_CONFIG_PATH,
# and MXLAPACK_INT.

# OpenBLAS provides lapack. On Ubuntu, the packages do not provlide lapack.h, so you need to install lapacke or lapacke64
# but still link only against openblas.
#
# On Ubuntu 22.04 installing openblas64 with a command like:
#
# $ sudo apt install libopenblas64-dev libopenblas64-serial-dev libopenblas64-openmp-dev
#
# installs the version(s) of openblas in threading-specific directories
# and we must set the pkg config path based on the choice.
#
# On systems with just openblas`, this will be ignored.
#
ifeq ($(findstring openblas, $(USE_BLAS_FROM)),openblas)

    ifeq ($(USE_BLAS_FROM),openblas64-serial)
       OPENBLAS_BLAS_PKG_CONFIG_NAME ?= openblas
       OPENBLAS_LAPACK_PKG_CONFIG_NAME ?= lapack
       OPENBLAS_PKG_CONFIG_PATH ?= /usr/lib/x86_64-linux-gnu/openblas64-serial/pkgconfig
       MXLAPACK_INT ?= long
    else ifeq ($(USE_BLAS_FROM),openblas64-openmp)
       OPENBLAS_BLAS_PKG_CONFIG_NAME ?= openblas
       OPENBLAS_LAPACK_PKG_CONFIG_NAME ?= lapack
       OPENBLAS_PKG_CONFIG_PATH ?= /usr/lib/x86_64-linux-gnu/openblas64-openmp/pkgconfig
       MXLAPACK_INT ?= long
    else ifeq ($(USE_BLAS_FROM),openblas64-pthread)
       OPENBLAS_BLAS_PKG_CONFIG_NAME ?= openblas
       OPENBLAS_LAPACK_PKG_CONFIG_NAME ?= lapack
       OPENBLAS_PKG_CONFIG_PATH ?= /usr/lib/x86_64-linux-gnu/openblas64-pthread/pkgconfig
       MXLAPACK_INT ?= long
    else ifeq ($(USE_BLAS_FROM),openblas-serial)
       OPENBLAS_BLAS_PKG_CONFIG_NAME ?= openblas
       OPENBLAS_LAPACK_PKG_CONFIG_NAME ?= lapack
       OPENBLAS_PKG_CONFIG_PATH ?= /usr/lib/x86_64-linux-gnu/openblas-serial/pkgconfig
       MXLAPACK_INT ?= int
    else ifeq ($(USE_BLAS_FROM),openblas-openmp)
       OPENBLAS_BLAS_PKG_CONFIG_NAME ?= openblas
       OPENBLAS_LAPACK_PKG_CONFIG_NAME ?= lapack
       OPENBLAS_PKG_CONFIG_PATH ?= /usr/lib/x86_64-linux-gnu/openblas-openmp/pkgconfig
       MXLAPACK_INT ?= int
    else ifeq ($(USE_BLAS_FROM),openblas-pthread)
       OPENBLAS_BLAS_PKG_CONFIG_NAME ?= openblas
       OPENBLAS_LAPACK_PKG_CONFIG_NAME ?= lapack
       OPENBLAS_PKG_CONFIG_PATH ?= /usr/lib/x86_64-linux-gnu/openblas-pthread/pkgconfig
       MXLAPACK_INT ?= int
    else ifeq ($(USE_BLAS_FROM),openblas64)
       #use pkg-config defaults
       OPENBLAS_BLAS_PKG_CONFIG_NAME ?= openblas64
       OPENBLAS_LAPACK_PKG_CONFIG_NAME ?= lapack64
       MXLAPACK_INT ?= long
    else
       #use pkg-config defaults
       OPENBLAS_PKG_CONFIG_PATH ?= $(PKG_CONFIG_PATH)
       OPENBLAS_BLAS_PKG_CONFIG_NAME ?= openblas
       OPENBLAS_LAPACK_PKG_CONFIG_NAME ?= lapack
       MXLAPACK_INT ?= int
    endif

    BLAS_INCLUDES ?= -DMXLIB_OPENBLAS -DMXLAPACK_INT=$(MXLAPACK_INT) $(shell PKG_CONFIG_PATH="$(OPENBLAS_PKG_CONFIG_PATH)" pkg-config --cflags "$(OPENBLAS_BLAS_PKG_CONFIG_NAME)") $(shell PKG_CONFIG_PATH="$(OPENBLAS_PKG_CONFIG_PATH)" pkg-config --cflags "$(OPENBLAS_LAPACK_PKG_CONFIG_NAME)")
    BLAS_LDFLAGS ?= $(shell PKG_CONFIG_PATH="$(OPENBLAS_PKG_CONFIG_PATH)" pkg-config --libs-only-L "$(OPENBLAS_BLAS_PKG_CONFIG_NAME)")
    BLAS_LDLIBS ?= $(shell PKG_CONFIG_PATH="$(OPENBLAS_PKG_CONFIG_PATH)" pkg-config --libs-only-l "$(OPENBLAS_BLAS_PKG_CONFIG_NAME)")


    ifeq ($(strip $(BLAS_LDFLAGS)),)
       $(error OpenBLAS ($(USE_BLAS_FROM)) not found by pkg-config)
    endif

endif

ifeq ($(USE_BLAS_FROM),atlas)
    ATLAS_BLAS_PKG_CONFIG_NAME ?= blas-atlas
    ATLAS_LAPACK_PKG_CONFIG_NAME ?= lapack-atlas
    MXLAPACK_INT ?= int

    BLAS_INCLUDES ?= -DMXLIB_OPENBLAS -DMXLAPACK_INT=$(MXLAPACK_INT) $(shell pkg-config --cflags "$(ATLAS_BLAS_PKG_CONFIG_NAME)") $(shell pkg-config --cflags "$(ATLAS_LAPACK_PKG_CONFIG_NAME)")
    BLAS_LDFLAGS ?= $(shell pkg-config --libs-only-L "$(ATLAS_BLAS_PKG_CONFIG_NAME)") $(shell pkg-config --libs-only-L "$(ATLAS_LAPACK_PKG_CONFIG_NAME)")
    BLAS_LDLIBS ?= $(shell pkg-config --libs-only-l "$(ATLAS_BLAS_PKG_CONFIG_NAME)") $(shell pkg-config --libs-only-l "$(ATLAS_LAPACK_PKG_CONFIG_NAME)")

    ifeq ($(strip $(BLAS_LDFLAGS)),)
       $(error ATLAS ($(USE_BLAS_FROM)) not found by pkg-config)
    endif

endif

ifeq ($(USE_FFT_FROM),fftw)
    
    ifeq ($(UNAME),Darwin)
        #Order matters, _threads first.
        FFT_LDLIBS ?= -lfftw3_threads -lfftw3f_threads -lfftw3l_threads -lfftw3 -lfftw3f -lfftw3l    
    endif
    ifeq ($(UNAME),Linux)
        FFT_LDLIBS ?= -lfftw3 -lfftw3f  -lfftw3l
    endif
endif

FITS_LIB = -lcfitsio

BOOST_LIB = -lboost_system -lboost_filesystem

GSL_LIB = -lgsl

ifeq ($(NEED_MILK),yes)

	MILK_INCLUDES ?= -I/usr/local/milk/include
	MILK_LIB ?= -L/usr/local/milk/lib -lImageStreamIO

	INCLUDES += $(MILK_INCLUDES)

	CFLAGS += -DMXLIB_MILK
	CXXFLAGS += -DMXLIB_MILK
endif


EXTRA_LDLIBS ?= $(FITS_LIB) $(BOOST_LIB) $(GSL_LIB) $(SOFA_LIB) $(MILK_LIB)

ifneq ($(UNAME),Darwin)
    EXTRA_LDLIBS += -lrt
endif

EXTRA_LDFLAGS += -L$(PREFIX)/lib

#FFTW before BLAS so MKL doesn't hijack certain functions
#FFTW:
EXTRA_LDLIBS += $(FFT_LDLIBS)
EXTRA_LDFLAGS += $(FFT_LDFLAGS)

#BLAS:
INCLUDES += $(BLAS_INCLUDES)
EXTRA_LDLIBS += $(BLAS_LDLIBS)
EXTRA_LDFLAGS += $(BLAS_LDFLAGS)

# Eigen
# -I/path/to/folder containing "Eigen" directory
EIGEN_CFLAGS ?= $(shell pkg-config eigen3 --cflags)
INCLUDES += $(EIGEN_CFLAGS)

LDLIBS += $(EXTRA_LDLIBS)
LDFLAGS += $(EXTRA_LDFLAGS)

CFLAGS += $(INCLUDES) $(OPTIMIZE)
CXXFLAGS += $(INCLUDES) $(OPTIMIZE) 

#This is needed to force use of g++ for linking
LINK.o = $(LINK.cc)

ifeq ($(NEED_CUDA),yes)

   CUDA_INCLUDES ?= -I/usr/local/cuda/include/

   CXXFLAGS += -DEIGEN_NO_CUDA $(CUDA_INCLUDES) -DMXLIB_CUDA


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

   HOST_COMPILER ?= g++
   NVCC          := nvcc -ccbin $(HOST_COMPILER)

   # internal flags
   NVCCFLAGS   := -m${TARGET_SIZE}
   NVCCFLAGS   +=  -DEIGEN_NO_CUDA -DMXLIB_MKL
   NVCCFLAGS   +=  ${NVCCARCH}

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

   #Finally we define the cuda libs for linking
   CUDA_LIBS ?= -L/usr/local/cuda/lib64/ -lcudart -lcublas -lcufft -lcurand

else
   CUDA_LIBS ?=
endif
