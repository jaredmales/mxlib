-include ../mk/MxLib.mk
-include ../mk/Common.mk


OBJS = \
	mxlib.o \
	mxError.o \
	ao/analysis/clGainOpt.o \
	ao/analysis/aoAtmosphere.o \
	ao/analysis/aoSystem.o \
	ao/analysis/fourierTemporalPSD.o \
	ao/analysis/aoPSDs.o \
	app/application.o \
	app/appConfigurator.o \
	app/clOptions.o \
	app/ini.o \
	app/optionparser.o \
	ioutils/fileUtils.o \
	ioutils/stringUtils.o \
	ioutils/fits/fitsUtils.o \
	ioutils/fits/fitsHeaderCard.o \
	ioutils/fits/fitsHeader.o \
	ioutils/textTable.o \
	improc/ADIDerotator.o \
	improc/ADIobservation.o \
	improc/HCIobservation.o \
	improc/KLIPreduction.o \
	improc/imageUtils.o \
	ipc/processInterface.o \
	ipc/sharedMemSegment.o \
	math/templateBLAS.o \
	math/templateLapack.o \
	math/fit/templateLevmar.o \
	math/fit/fitGaussian.o \
	math/fft/fft.o \
	math/fft/fftwEnvironment.o \
	math/fft/fftwTemplates.o \
	math/func/bessel.o \
	math/func/factorial.o \
	math/func/gamma.o \
	math/func/jinc.o \
	math/func/legendre.o \
	math/func/moffat.o \
	math/func/precision.o \
	math/func/sign.o \
	math/plot/gnuPlot.o \
	sigproc/zernike.o \
	sys/environment.o \
	sys/timeUtils.o

ifeq ($(NEED_CUDA),yes)

CUDA_OBJS = \
	math/cuda/templateCublas.o \
	math/cuda/templateCufft.o \
	math/cuda/templateCurand.o
endif

VENDOROBJS = \
	vendor/levmar-2.6/lm.o \
	vendor/levmar-2.6/Axb.o \
	vendor/levmar-2.6/misc.o \
	vendor/levmar-2.6/lmlec.o \
	vendor/levmar-2.6/lmbc.o \
	vendor/levmar-2.6/lmblec.o \
	vendor/levmar-2.6/lmbleic.o

CURRENT_SOFA = vendor/sofa/20210125

all: libmxlib.$(LIB_SUFFIX)

ao/analysis/aoSystem.o: ../include/ao/analysis/aoSystem.hpp
ao/analysis/fourierTemporalPSD.o: ../include/ao/analysis/fourierTemporalPSD.hpp
improc/HCIobservation.o: ../include/improc/HCIobservation.hpp
ioutils/fileUtils.o: ../include/ioutils/fileUtils.hpp
math/cuda/templateCublas.o: ../include/math/cuda/templateCublas.hpp
math/cuda/templateCufft.o: ../include/math/cuda/templateCufft.hpp
math/cuda/templateCurand.o: ../include/math/cuda/templateCurand.hpp
math/func/bessel.o: ../include/math/func/bessel.hpp
math/func/legendre.o: ../include/math/func/legendre.hpp
app/application.o: ../include/app/application.hpp ../include/app/appConfigurator.hpp ../include/app/clOptions.hpp
app/appConfigurator.o: ../include/app/appConfigurator.hpp ../include/app/clOptions.hpp
app/clOptions.o: ../include/app/clOptions.hpp
app/ini.o: ../include/app/ini.hpp
app/optionparser.o: ../include/app/optionparser/optionparser.h

#mxlib_uncomp_version

libmxlib.$(LIB_SUFFIX): sofa $(OBJS) $(CUDA_OBJS)
	$(MAKE) -C vendor/levmar-2.6/
ifeq ($(UNAME),Darwin)
	$(LINK.o) -dynamiclib -install_name @rpath/libmxlib.dylib -rpath @loader_path $(VENDOROBJS) $(OBJS) $(CUDA_OBJS) -o libmxlib.dylib $(LDFLAGS) $(LDLIBS) $(CUDA_LIBS)
else
	$(LINK.o) -shared $(VENDOROBJS) $(OBJS) $(CUDA_OBJS) -o $@ $(LDFLAGS) $(LDLIBS) $(CUDA_LIBS)
endif

sofa:
	$(MAKE) -C $(CURRENT_SOFA)/c/src/ install

install: libmxlib.$(LIB_SUFFIX)
	install -d $(LIB_PATH)
	install libmxlib.$(LIB_SUFFIX) $(LIB_PATH)

.PHONY: clean
clean:
	rm -f $(OBJS) $(CUDA_OBJS)
	$(MAKE) -C vendor/levmar-2.6/ clean
	$(MAKE) -C $(CURRENT_SOFA)/c/src/ clean

.PHONY: realclean
realclean: clean
	rm -f libmxlib.$(LIB_SUFFIX)
	rm -f $(OBJS)
	$(MAKE) -C vendor/levmar-2.6/ cleanall

# .PHONY: printvars
# printvars:
#     @$(foreach V,$(sort $(.VARIABLES)),$(if $(filter-out environment% default automatic,$(origin $V)),$(warning $V=$($V) ($(value $V)))))
