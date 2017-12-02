# Common path and make variable definitions
#
# NOTE: This file should only be edited in mxlib/local, not in the root mxlib directory.
#
-include local/Common.mk

UNAME ?= $(shell uname)
ifeq ($(UNAME),Darwin)
	CFLAGS += -D_BSD_SOURCE
	CXXFLAGS += -D_BSD_SOURCE
else
	CFLAGS += -D_OPEN_SOURCE=700
	CXXFLAGS += -D_OPEN_SOURCE=700
endif
PREFIX ?= $(HOME)
BIN_PATH ?= $(PREFIX)/bin
LIB_PATH ?= $(PREFIX)/lib
INCLUDE_PATH ?= $(PREFIX)/include
LIB_SOFA ?= $(LIB_PATH)/libsofa_c.a
ARFLAGS ?= rvs

INCLUDES += -I$(INCLUDE_PATH)