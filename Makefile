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

OBJS = kepler.o \
       astrodyn.o \
       msgq.o \
       mxlib.o\
       sharedmem_segment.o \
       sharedMemSegment.o \
       process_interface.o \
       ds9_interface.o \


INC_TO_INSTALL = ADIobservation.hpp \
                 airy.hpp \
                 app \
                 appConfigurator.hpp \
                 application.hpp \
                 astroconstants.h \
                 astrodyn.hpp \
                 astroFilter.hpp \
                 astroSpectrum.hpp \
                 astrotypes.h \
                 autocorrelation.hpp \
                 binVector.hpp \
                 ds9_interface.h \
                 eigenImage.hpp \
                 eigenCube.hpp \
                 eigenUtils.hpp \
                 environment.hpp \
                 fft.hpp \
                 fftwTemplates.hpp \
                 fileUtils.hpp \
                 fitsUtils.hpp \
                 fitsFile.hpp \
                 fitsHeader.hpp \
                 fitsHeaderCard.hpp \
                 fourierModes.hpp \
                 fraunhoferImager.hpp \
                 gaussian.hpp \
                 geo.h \
                 gnuPlot.hpp \
                 gramSchmidt.hpp \
                 gslInterpolation.hpp \
       		 HCIobservation.hpp \
                 IPC.h \
                 imageFilters.hpp \
       		 imageMasks.hpp \
       		 imagePads.hpp \
                 imageTransforms.hpp \
                 imagingArray.hpp \
                 imagingUtils.hpp \
                 jinc.hpp \
                 KLIPreduction.hpp \
                 kepler.hpp \
                 levmarInterface.hpp \
                 lyotCoronagraph.hpp \
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
                 sharedMemSegment \
                 sigmoid.hpp \
                 signalWindows.hpp \
                 stringUtils.hpp \
                 templateBLAS.hpp \
                 templateLapack.hpp \
                 templateLevmar.hpp \
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
sharedMemSegment.o: include/sharedMemSegment
process_interface.o: include/process_interface.h
ds9_interface.o: include/ds9_interface.h
astrodyn.o: include/astrodyn.hpp
kepler.o: include/kepler.hpp

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
