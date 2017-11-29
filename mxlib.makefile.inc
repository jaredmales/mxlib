#Definitions for building mxlib
#
# NOTE: This file should only be edited in mxlib/local, not in the root mxlib directory.
#
INSTALL_PATH = $(HOME)
INCLUDE_PATH = $(INSTALL_PATH)/include/mx
LIB_PATH = $(INSTALL_PATH)/lib
BIN_PATH = $(INSTALL_PATH)/bin 

LIB_SOFA = $(LIB_PATH)/libsofa_c.a

OPTIMIZE = -O3

CXX = g++
AR = ar -r
RANLIB = ar -s

#must include path to the include directory, and to sofa
INCLUDE = -Iinclude -I$(HOME)/include

CFLAGS += --std=c99 -D_XOPEN_SOURCE=600  -fPIC
CXXFLAGS += --std=c++14 -D_XOPEN_SOURCE=600 -fPIC

