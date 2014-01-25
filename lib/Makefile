INCLUDE_PATH = $(HOME)/include
LIB_PATH = $(HOME)/lib

OPTIMIZE = -O3

CPP = g++
AR = ar -r
RANLIB = ar -s

CPPFLAGS += -std=c++0x -ggdb
.c.o:
	$(CC) $(OPTIMIZE) $(CFLAGS) -Iinclude -c $<

.cpp.o:
	$(CPP) $(OPTIMIZE) $(CPPFLAGS) -Iinclude -c $<

# programs to be made
TARGETS = libcda
OBJS = readcolumns.o kepler.o astrodyn.o

INC_TO_INSTALL = astroconstants.h astrodyn.h astrotypes.h kepler.h geo.h readcolumns.h
VMOP_TO_INSTALL = MatrixOp MMatrix MMatrix1 MMatrix2 MMatrixSlice MMatrixView1 MMatrixView2 vmopFitsUtils

all: $(TARGETS) 


libcda: $(OBJS) vmop/MatrixOp.o vmop/vmopFitsUtils.o sofa/libsofa_c.a
	$(AR) libcda.a $(OBJS) MatrixOp.o vmopFitsUtils.o sofa/*.o
	$(RANLIB) libcda.a

sofa/libsofa_c.a:
	$(MAKE) -C ./sofa all

install:
	install -d $(INCLUDE_PATH)
	install libcda.a $(LIB_PATH)
	for file in ${INC_TO_INSTALL}; do \
	 (install include/$$file $(INCLUDE_PATH)/cDA) || break; \
	done
	for file in ${VMOP_TO_INSTALL}; do \
	 (install vmop/$$file $(INCLUDE_PATH)/cDA) || break; \
	done
	install sofa/*.h $(INCLUDE_PATH)/cDA

clean:
	rm -f *.o *~
	rm -f libcda.a
	$(MAKE) -C ./sofa clean
