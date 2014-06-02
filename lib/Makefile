INSTALL_PATH = $(HOME)
INCLUDE_PATH = $(INSTALL_PATH)/include/mx
LIB_PATH = $(INSTALL_PATH)/lib

OPTIMIZE = -ggdb

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
       msgq.o \
       sharedmem_segment.o \
       sharedMemSegment.o \
       process_interface.o \
       ds9_interface.o \
       gnuplot_interface.o


INC_TO_INSTALL = astroconstants.h \
                 astrodyn.h \
                 astrotypes.h \
                 kepler.h \
                 geo.h \
                 readcolumns.h \
                 stringUtils \
                 timeUtils \
                 pout \
                 IPC.h \
                 msgq.h \
                 msgQ \
                 sharedmem_segment.h \
                 sharedMemSegment \
                 process_interface.h \
                 ds9_interface.h \
                 gnuplot_interface.h \
                 fitsUtils \
                 fitsFile \
                 fitsHeader \
                 fitsHeaderCard \
                 mxException
                 
VMOP_TO_INSTALL = MMatrix \
                  MMatrix1 \
                  MMatrix2 \
                  MMatrixSlice \
                  MMatrixView1 \
                  MMatrixView2

all: $(TARGETS) 

#dependencies:
msgq.o: include/IPC.h include/msgq.h
sharedmem_segment.o: include/IPC.h include/sharedmem_segment.h
sharedMemSegment.o: include/sharedMemSegment
process_interface.o: include/process_interface.h
ds9_interface.o: include/ds9_interface.h
gnuplot_interface.o: include/gnuplot_interface.h

libmxlib: $(OBJS)
	$(AR) libmxlib.a $(OBJS)
	$(RANLIB) libmxlib.a
	gcc -shared -Wl,-soname,libmxlib.so -o libmxlib.so $(OBJS) /home/jaredmales/lib/libsofa_c.a

install: libmxlib
	install -d $(INCLUDE_PATH)
	install -d $(LIB_PATH)
	install libmxlib.a $(LIB_PATH)
	install libmxlib.so $(LIB_PATH)
	for file in ${INC_TO_INSTALL}; do \
	 (install include/$$file $(INCLUDE_PATH)) || break; \
	done
	for file in ${VMOP_TO_INSTALL}; do \
	 (install vmop/$$file $(INCLUDE_PATH)) || break; \
	done

clean:
	rm -f *.o *~
	rm -f libmxlib.a
