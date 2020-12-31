-include ../mk/MxLib.mk
-include ../mk/Common.mk


OBJS = mxError.o \
       ao/analysis/clGainOpt.o \
       ioutils/fileUtils.o \
       ioutils/stringUtils.o \
       ioutils/fits/fitsUtils.o \
       ioutils/fits/fitsHeaderCard.o \
       ioutils/fits/fitsHeader.o \
       math/cuda/templateCublas.o \
       math/cuda/templateCufft.o \
       math/cuda/templateCurand.o \
       math/templateBLAS.o \
       math/templateLapack.o \
       math/fit/templateLevmar.o \
       math/fft/fft.o \
       math/fft/fftwEnvironment.o \
       math/fft/fftwTemplates.o \
       math/func/bessel.o \
       math/func/factorial.o \
       math/func/gamma.o \
       math/func/precision.o \
       math/func/sign.o \
       math/plot/gnuPlot.o \
       sys/environment.o \
       sys/timeUtils.o 
       

all: libmxlib.so

#Dependencies:
#mxError.o: ../include/mxError.hpp 
#ao/analysis/clGainOpt.o:   ../include/ao/analysis/clGainOpt.hpp 
#ioutils/fileUtils.o: ../include/ioutils/fileUtils.hpp 
#ioutils/stringUtils.o: ../include/ioutils/stringUtils.hpp 
#ioutils/fits/fitsUtils.o: ../include/ioutils/fits/fitsUtils.hpp      
#ioutils/fits/fitsHeaderCard.o: ../include/ioutils/fits/fitsHeaderCard.hpp
#ioutils/fits/fitsHeader.o: ../include/ioutils/fits/fitsHeader.hpp
#math/cuda/templateCublas.o: ../include/math/cuda/templateCublas.hpp
#math/cuda/templateCufft.o: ../include/math/cuda/templateCufft.hpp
#math/cuda/templateCurand.o: ../include/math/cuda/templateCurand.hpp
#math/templateBLAS.o: ../include/math/templateBLAS.hpp
#math/templateLapack.o: ../include/math/templateLapack.hpp
#math/fit/templateLevmar.o: ../include/math/fit/templateLevmar.hpp
#math/fft/fft.o: ../include/math/fft/fft.hpp
#math/fft/fftwEnvironment.o: ../include/math/fft/fftwEnvironment.hpp
#math/fft/fftwTemplates.o: ../include/math/fft/fftwTemplates.hpp
#math/func/bessel.o: ../include/math/func/bessel.hpp
#math/func/factorial.o: ../include/math/func/factorial.hpp
#math/func/gamma.o: ../include/math/func/gamma.hpp
#math/func/precision.o: ../include/math/func/precision.hpp
#math/func/sign.o: ../include/math/func/sign.hpp
#math/plot/gnuPlot.o: ../include/math/plot/gnuPlot.hpp
#sys/environment.o: ../include/sys/environment.hpp
#sys/timeUtils.o: ../include/sys/timeUtils.hpp
			

#mxlib_uncomp_version

libmxlib.so: $(OBJS) 
	$(CXX) -shared $(OBJS) -o $@

install: libmxlib.so
	install -d $(LIB_PATH)
	install libmxlib.so $(LIB_PATH)

.PHONY: clean
clean:
	rm -f libmxlib.so
	rm -f $(OBJS)
	rm -f *~
	
