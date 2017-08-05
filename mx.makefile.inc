#Definitions for building mxlib based applications
#
# NOTE: This file should only be edited in mxlib/local, not in the root mxlib directory.
#


INSTALL_DIR = $(HOME)
INCLUDE_PATH = $(INSTALL_DIR)/include
LIB_PATH = $(INSTALL_DIR)/lib
BIN_PATH = $(INSTALL_DIR)/bin 

FFT_LIB = 
BLASPACK_INC = -DMXLIB_MKL_ONLY -DMKL_ILP64 -m64 -I${MKLROOT}/include
BLASPACK_LIB = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl

#FFT_LIB = -lfftw3 -lfftw3f
#BLASPACK_INC = -I/usr/local/atlas/include
#BLASPACK_LIB = -L/usr/local/atlas/lib/ -llapack -lf77blas -lcblas -latlas -lgfortran

MXLIB_EXLIBS = -lsofa_c -llevmar -L/usr/lib64/ -lcfitsio -lrt -L/usr/local/lib -lboost_system -lboost_filesystem $(FFT_LIB) -lgsl $(BLASPACK_LIB) 

OPTIMIZE = -O3 -fopenmp -ffast-math

CPP = g++

CFLAGS += --std=c99 -D_XOPEN_SOURCE=600  -fPIC
CPPFLAGS += --std=c++14 -D_XOPEN_SOURCE=600 -fPIC

INCLUDE = -I$(INCLUDE_PATH) $(BLASPACK_INC)

.c.o:
	$(CC) $(OPTIMIZE) $(CFLAGS) $(INCLUDE) -c $<

.cpp.o:
	$(CPP) $(OPTIMIZE) $(CPPFLAGS) $(INCLUDE) -c $<


	
	
