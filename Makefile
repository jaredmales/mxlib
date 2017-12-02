##########################################
##                                      ##
##        Makefile for mxlib            ##
##                                      ##
##########################################

# Instruction:

# Step 1: run make -f Makefile.setup 
#   -- this creates the ./local directory, and copies mxlib.makefile.inc there.
#  
# Step 2: edit local/mxlib.makefile.inc 
# NOTE: in general, you should only edit ./local/mxlib.makefile.inc, not
# any files in mxlib/.  This is to preserve the git repository state.


include local/mxlib.makefile.inc 
UNAME := $(shell uname)
LIBNAME = libmxlib
ifeq ($(UNAME),Darwin)  # macOS
	LIBEXT = dylib
	LIBFLAGS = -dynamiclib -install_name "$(LIBNAME).$(LIBEXT)" -o $(LIBNAME).$(LIBEXT)
else
	LIBEXT = so
    LIBFLAGS = -Wl,-soname,$(LIBNAME).$(LIBEXT) -o $(LIBNAME).$(LIBEXT) -lrt
endif

.c.o:
	$(CC) $(OPTIMIZE) $(CFLAGS) $(INCLUDE) -c $<

.cpp.o:
	$(CXX) $(OPTIMIZE) $(CXXFLAGS) $(INCLUDE) -c $<

# programs to be made
TARGETS = libmxlib

OBJS = msgq.o \
       mxlib.o\
       sharedmem_segment.o \
       process_interface.o 

#astrodyn.o \

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

#dependencies:
msgq.o: include/IPC.h include/msgq.h
mxlib.o: include/mxlib.h include/mxlib_comp_version.h
sharedmem_segment.o: include/IPC.h include/sharedmem_segment.h
process_interface.o: include/process_interface.h
astrodyn.o: include/astrodyn.hpp

.PHONY: mxlib_comp_version
mxlib_comp_version:
	@sh ./gengithead.sh ./ ./include/mxlib_comp_version.h MXLIB_COMP

.PHONY: mxlib_uncomp_version
mxlib_uncomp_version:
	@sh ./gengithead.sh ./ ./include/mxlib_uncomp_version.h MXLIB_UNCOMP
	
libmxlib: mxlib_comp_version mxlib_uncomp_version $(OBJS) 
	$(AR) libmxlib.a $(OBJS)
	$(RANLIB) libmxlib.a 
	$(CC) -shared $(LIBFLAGS) $(OBJS) $(LIB_SOFA) -lc -rdynamic
	
install: libmxlib
	install -d $(INCLUDE_PATH)
	install -d $(LIB_PATH)
	install $(LIBNAME).a $(LIB_PATH)
	install $(LIBNAME).$(LIBEXT) $(LIB_PATH)
	install gengithead.sh $(BIN_PATH)
	for file in ${INC_TO_INSTALL}; do \
	 (cp -r include/$$file $(INCLUDE_PATH)) || break; \
	done
# 	for file in ${VMOP_TO_INSTALL}; do \
# 	 (install vmop/$$file $(INCLUDE_PATH)) || break; \
# 	done

clean:
	rm -f *.o *~
	rm -f libmxlib.a
	rm -f libmxlib.so
