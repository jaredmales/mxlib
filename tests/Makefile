-include ../local/MxLib.mk

-include ../mk/Common.mk

INCLUDES += -I../include

OBJS = testsMain.o \
       include/astro/astroDynamics_test.o \
       include/ioutils/fileUtils_test.o \
       include/math/cuda/templateCublas_test.o \
       include/math/templateBLAS_test.o \
       include/math/templateLapack_test.o \
       include/math/randomT_test.o \
       include/sigproc/psdUtils_test.o \
       include/sigproc/psdFilter_test.o \
       include/sigproc/zernike_test.o \
       include/sys/timeUtils_test.o \
       
       
       

all: mxlibTest

			

#mxlib_uncomp_version

mxlibTest: $(OBJS) 
	$(LINK.o)  -o $@ $(OBJS) ../source/libmxlib.so $(LDFLAGS) $(LDLIBS) $(CUDA_LIBS)


.PHONY: clean
clean:
	rm -f mxlibTest
	rm -f $(OBJS)
	rm -f *~
	
