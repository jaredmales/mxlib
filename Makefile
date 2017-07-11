INSTALL_PATH = $(HOME)
INCLUDE_PATH = $(INSTALL_PATH)/include/mx
LIB_PATH = $(INSTALL_PATH)/lib
BIN_PATH = $(HOME)/bin 

LIB_SOFA = $(LIB_PATH)/libsofa_c.a

OPTIMIZE = -O3

CPP = g++
AR = ar -r
RANLIB = ar -s

#must include path to the include directory, and to sofa
INCLUDE = -Iinclude -I$(HOME)/include

CFLAGS += --std=c99 -D_XOPEN_SOURCE=600  -fPIC
CPPFLAGS += --std=c++0x -D_XOPEN_SOURCE=600 -fPIC

.c.o:
	$(CC) $(OPTIMIZE) $(CFLAGS) $(INCLUDE) -c $<

.cpp.o:
	$(CPP) $(OPTIMIZE) $(CPPFLAGS) $(INCLUDE) -c $<

# programs to be made
TARGETS = libmxlib

OBJS = astrodyn.o \
       msgq.o \
       mxlib.o\
       sharedmem_segment.o \
       process_interface.o 


INC_TO_INSTALL = ao \
                 app \
                 astro \
                 astrodyn.hpp \
                 astroFilter.hpp \
                 astroSpectrum.hpp \
                 astrotypes.h \
                 autocorrelation.hpp \
                 binVector.hpp \
                 eigenUtils.hpp \
                 environment.hpp \
                 fft.hpp \
                 fftwTemplates.hpp \
                 fileUtils.hpp \
                 fitGaussian.hpp \
                 fourierModes.hpp \
                 gaussian.hpp \
                 geo.hpp \
                 gnuPlot.hpp \
                 gramSchmidt.hpp \
                 gslInterpolation.hpp \
                 IPC.h \
                 improc \
                 wfp \
                 imagingArray.hpp \
                 jinc.hpp \
                 levmarInterface.hpp \
                 logistic.hpp \
                 msgq.h \
                 msgQ \
                 mxException.hpp \
                 mxError.hpp \
                 mxlib.h\
                 mxlib_uncomp_version.h\
                 ompLoopWatcher.hpp \
                 pout.hpp \
                 process_interface.h \
                 psdFilter.hpp \
                 psdUtils.hpp \
                 readColumns.hpp \
                 roots.hpp \
                 sharedmem_segment.h \
                 signalWindows.hpp \
                 stringUtils.hpp \
                 templateBLAS.hpp \
                 templateLapack.hpp \
                 templateLevmar.hpp \
                 textTable.hpp \
                 trueFalseT.hpp \
                 timeUtils.hpp \
		 randomSeed.hpp \
		 randomT.hpp \
		 vectorUtils.hpp  \
		 zernike.hpp

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
	gcc -shared -Wl,-soname,libmxlib.so -o libmxlib.so $(OBJS) $(LIB_SOFA) -lrt -lc -rdynamic
	
install: libmxlib
	install -d $(INCLUDE_PATH)
	install -d $(LIB_PATH)
	install libmxlib.a $(LIB_PATH)
	install libmxlib.so $(LIB_PATH)
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
