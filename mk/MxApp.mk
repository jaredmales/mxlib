# Makefile for building mxlib based applications
#
SELF_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
-include $(SELF_DIR)/../local/MxApp.mk
include $(SELF_DIR)/../mk/Common.mk


# Single-file app name can be supplied as `TARGET=`,
# or `t=` for short
TARGET ?= $(t)

#allow passing of .o, .hpp, or any extension, or no extension
BASENAME=$(basename $(TARGET))
PATHNAME=$(dir $(TARGET))
OBJNAME = $(PATHNAME)$(notdir $(BASENAME)).o
TARGETNAME = $(PATHNAME)$(notdir $(BASENAME))

#Check if we want a target-specific git version header
ifeq ($(GIT_VERSION),yes)
   #change target name to uppercase
   GIT_VERSION_DEF = $(shell echo $(TARGET) | tr a-z A-Z)_GIT
   PRE_TARGETS += git_version
   GIT_VERSION_FILE = $(TARGETNAME)_git_version.h
   GIT_PATH ?= ./
endif



all: tshoot $(PRE_TARGETS) $(TARGETNAME) $(OTHER_OBJS)

$(TARGETNAME): $(OBJNAME) $(OTHER_OBJS)
	$(LINK.o)  -o $(TARGETNAME) $(OBJNAME) $(OTHER_OBJS) -lmxlib $(LDFLAGS) $(LDLIBS) $(CUDA_LIBS) $(OTHER_LIBS)


install: all
	install -d $(BIN_PATH)
	install $(TARGET) $(BIN_PATH)

.PHONY: tshoot
tshoot:
	@echo mxlib: compiling $(BASENAME)

.PHONY: git_version
git_version:
	@bash gengithead.sh $(GIT_PATH) $(TARGETNAME)_git_version.h $(GIT_VERSION_DEF)
	
.PHONY: clean
clean:
	rm -f $(TARGETNAME) $(OBJNAME) $(GIT_VERSION_FILE)
	rm -f *.o 
	rm -f *~
