#Definitions for building mxlib based applications
#
# NOTE: This file should only be edited in mxlib/local, not in the root mxlib directory.
#


INSTALL_DIR = $(HOME)
INCLUDE_PATH = $(INSTALL_DIR)/include
LIB_PATH = $(INSTALL_DIR)/lib
BIN_PATH = $(INSTALL_DIR)/bin 

MXLIB_EXLIBS = -lsofa_c -llevmar -L/usr/lib64/ -lcfitsio -lrt -L/usr/local/lib -lboost_system -lboost_filesystem -L/usr/local/atlas/lib/ -ltatlas -lfftw3 -lfftw3f -lgsl -lgslcblas

OPTIMIZE = -O3 -fopenmp -ffast-math

CPP = g++

CFLAGS += --std=c99 -D_XOPEN_SOURCE=600  -fPIC
CPPFLAGS += --std=c++14 -D_XOPEN_SOURCE=600 -fPIC

INCLUDE = -I../include -I$(INCLUDE_PATH) -I/usr/local/atlas/include

.c.o:
	$(CC) $(OPTIMIZE) $(CFLAGS) $(INCLUDE) -c $<

.cpp.o:
	$(CPP) $(OPTIMIZE) $(CPPFLAGS) $(INCLUDE) -c $<


	
	
