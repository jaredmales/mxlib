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

OBJS = readcolumns.o \
       kepler.o \
       astrodyn.o \
       fileUtils.o \
       msgq.o \
       mxlib.o\
       sharedmem_segment.o \
       sharedMemSegment.o \
       process_interface.o \
       ds9_interface.o \
       gnuplot_interface.o \


INC_TO_INSTALL = app \
                 appConfigurator.hpp \
                 application.hpp \
                 astroconstants.h \
                 astrodyn.hpp \
                 astrotypes.h \
                 kepler.hpp \
                 geo.h \
                 readcolumns.h \
                 stringUtils.hpp \
                 fileUtils.hpp \
                 gslInterpolation.hpp \
                 pout.hpp \
                 IPC.h \
                 msgq.h \
                 msgQ \
                 mxlib.h\
                 mxlib_uncomp_version.h\
                 sharedmem_segment.h \
                 sharedMemSegment \
                 process_interface.h \
                 ds9_interface.h \
                 gnuplot_interface.h \
                 fileUtils.hpp \
                 fitsUtils.hpp \
                 fitsFile.hpp \
                 fitsHeader.hpp \
                 fitsHeaderCard.hpp \
                 eigenImage.hpp \
                 eigenCube.hpp \
                 eigenUtils.hpp \
                 imageFilters.hpp \
                 imageTransforms.hpp \
                 templateBLAS.hpp \
                 templateLapack.hpp \
                 templateLevmar.hpp \
                 timeUtils.hpp \
                 mxException.hpp \
		 HCIobservation.hpp \
		 ADIobservation.hpp \
	         KLIPreduction.hpp \
		 imageMasks.hpp
                 
# VMOP_TO_INSTALL = MMatrix \
#                   MMatrix1 \
#                   MMatrix2 \
#                   MMatrixSlice \
#                   MMatrixView1 \
#                   MMatrixView2

all: $(TARGETS) 

#dependencies:
fileUtils.o: include/fileUtils.hpp
msgq.o: include/IPC.h include/msgq.h
mxlib.o: include/mxlib.h include/mxlib_comp_version.h
sharedmem_segment.o: include/IPC.h include/sharedmem_segment.h
sharedMemSegment.o: include/sharedMemSegment
process_interface.o: include/process_interface.h
ds9_interface.o: include/ds9_interface.h
gnuplot_interface.o: include/gnuplot_interface.h
astrodyn.o: include/astrodyn.hpp
kepler.o: include/kepler.hpp

.PHONY: mxlib_comp_version
mxlib_comp_version:
	@sh ./gengithead.sh ../ ./include/mxlib_comp_version.h MXLIB_COMP

.PHONY: mxlib_uncomp_version
mxlib_uncomp_version:
	@sh ./gengithead.sh ../ ./include/mxlib_uncomp_version.h MXLIB_UNCOMP
	
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
