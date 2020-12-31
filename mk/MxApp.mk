# Definitions for building mxlib based applications
#
SELF_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
-include $(SELF_DIR)/../local/MxApp.mk
include $(SELF_DIR)/../mk/Common.mk


# Single-file app name can be supplied as `TARGET=`,
# or `t=` for short
TARGET ?= $(t)

#Check if we want a target specific git version header
ifeq ($(GIT_VERSION),yes)
   #change target name to uppercase
   GIT_VERSION_DEF = $(shell echo $(TARGET) | tr a-z A-Z)_GIT
   PRE_TARGETS += git_version
endif


all: $(PRE_TARGETS) $(TARGET) $(OTHER_OBJS) 

$(TARGET):  $(TARGET).o  $(OTHER_OBJS)
	$(LINK.o)  -o $(TARGET) $(TARGET).o $(OTHER_OBJS) -lmxlib $(LDFLAGS) $(LDLIBS) $(CUDA_LIBS)

#endif

install: all
	install -d $(BIN_PATH)
	install $(TARGET) $(BIN_PATH)

.PHONY: git_version
git_version:
	@bash gengithead.sh ./ ./$(TARGET)_git_version.h $(GIT_VERSION_DEF)
	
.PHONY: clean
clean:
	rm -f $(TARGET) $(TARGET)_git_version.h
	rm -f *.o 
	rm -f *~
