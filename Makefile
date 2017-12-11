##########################################
##                                      ##
##          Makefile for mxlib          ##
##                                      ##
##########################################

# Are you here to customize the build?
# 
# Option 1: Override the options on the command line
#     e.g. `make install PREFIX=$HOME`
# Option 2: `make setup`
#     This creates ./local/{Common,MxLib,MxApp}.mk and lists options.
#     Define those options in the ./local Makefiles. These files are
#     *not* tracked in git, so your tweaks won't cause warnings about
#     uncommitted changes.

include mk/Common.mk
include mk/MxLib.mk

ifeq ($(UNAME),Darwin)  # macOS
	LIBEXT = dylib
	SHAREDLIBFLAGS = -dynamiclib -install_name "libmxlib.$(LIBEXT)" -o libmxlib.$(LIBEXT)
else
	LIBEXT = so
	SHAREDLIBFLAGS = -Wl,-soname,libmxlib.$(LIBEXT) -o libmxlib.$(LIBEXT) -lrt
endif

CFLAGS += -std=c99 -fPIC $(OPTIMIZE) $(INCLUDES)
CXXFLAGS += -std=c++14 -fPIC $(OPTIMIZE) $(INCLUDES)

# Programs to be made:
TARGETS = libmxlib

OBJS = msgq.o \
       mxlib.o\
       sharedmem_segment.o 
#\
#process_interface.o 

INC_TO_INSTALL = ao \
                 app \
                 astro \
                 fft \
                 improc \
                 math \
                 meta \
                 sigproc \
                 wfp \
                 astrodyn.hpp \
                 astrotypes.h \
                 binVector.hpp \
                 eigenUtils.hpp \
                 environment.hpp \
                 fileUtils.hpp \
                 gnuPlot.hpp \
                 gslInterpolation.hpp \
                 IPC.h \
                 imagingArray.hpp \
                 msgq.h \
                 msgQ \
                 mxException.hpp \
                 mxError.hpp \
                 mxlib.h\
                 mxlib_uncomp_version.h\
                 ompLoopWatcher.hpp \
                 pout.hpp \
                 process_interface.h \
                 readColumns.hpp \
                 sharedmem_segment.h \
                 stringUtils.hpp \
                 textTable.hpp \
                 timeUtils.hpp \
				 randomSeed.hpp \
				 randomT.hpp 

all: $(TARGETS) 

# Dependencies:
msgq.o: include/IPC.h include/msgq.h
mxlib.o: include/mxlib.h include/mxlib_comp_version.h
sharedmem_segment.o: include/IPC.h include/sharedmem_segment.h
#process_interface.o: include/process_interface.h
#astrodyn.o: include/astrodyn.hpp

.PHONY: mxlib_comp_version
mxlib_comp_version:
	@sh ./gengithead.sh ./ ./include/mxlib_comp_version.h MXLIB_COMP

.PHONY: mxlib_uncomp_version
mxlib_uncomp_version:
	@sh ./gengithead.sh ./ ./include/mxlib_uncomp_version.h MXLIB_UNCOMP

.PHONY: setup
setup:
	@for file in ./local/*.example.mk; do \
		dest=$$(echo $$file | sed 's/.example//'); \
		if [ ! -e $$dest ]; then cp -v $$file $$dest; fi \
	done
	@echo "***\nBuild settings available in local/Common.mk\n***"
	@grep "?=" mk/Common.mk || true
	@echo "***"
	@echo "Build settings available in local/MxLib.mk\n***"
	@grep "?=" mk/MxLib.mk || true
	@echo "***"
	@echo "Build settings available in local/MxApp.mk\n***"
	@grep  "?=" mk/MxApp.mk || true
	@echo "***"

libmxlib: mxlib_comp_version mxlib_uncomp_version $(OBJS) 
	$(AR) $(ARFLAGS) libmxlib.a $(OBJS)
	$(CC) -shared $(SHAREDLIBFLAGS) $(OBJS) $(LIB_SOFA) -lc -rdynamic
	
install: libmxlib
	install -d $(INCLUDE_PATH)/mx
	install -d $(LIB_PATH)
	install libmxlib.a $(LIB_PATH)
	install libmxlib.$(LIBEXT) $(LIB_PATH)
	install gengithead.sh $(BIN_PATH)
	for file in ${INC_TO_INSTALL}; do \
	  (cp -r include/$$file $(INCLUDE_PATH)/mx) || break; \
	done

clean:
	rm -f *.o *~
	rm -f libmxlib.a
	rm -f libmxlib.$(LIBEXT)
