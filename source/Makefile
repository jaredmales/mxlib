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
	sys/gitRepo.o \
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

all: libmxlib.$(LIB_SUFFIX)

mxlib.o: .FORCE
mxError.o: ../include/mxError.hpp ../include/mxException.hpp
ao/analysis/aoAtmosphere.o: ../include/ao/analysis/aoAtmosphere.hpp 
ao/analysis/aoPSDs.o: ../include/ao/analysis/aoPSDs.hpp
ao/analysis/aoSystem.o: ../include/ao/analysis/aoSystem.hpp ao/analysis/aoAtmosphere.o ao/analysis/aoPSDs.o ../include/ao/analysis/aoWFS.hpp
ao/analysis/clGainOpt.o: ../include/ao/analysis/clGainOpt.hpp
ao/analysis/clAOLinearPredictor.o: ../include/ao/analysis/clAOLinearPredictor.hpp ao/analysis/clGainOpt.o
ao/analysis/fourierTemporalPSD.o: ../include/ao/analysis/fourierTemporalPSD.hpp ../include/math/gslInterpolation.hpp ../include/sigproc/psdVarMean.hpp ao/analysis/aoSystem.o ao/analysis/clGainOpt.o ao/analysis/clAOLinearPredictor.o
app/application.o: ../include/app/application.hpp ../include/app/appConfigurator.hpp ../include/app/clOptions.hpp
app/appConfigurator.o: ../include/app/appConfigurator.hpp ../include/app/clOptions.hpp
app/clOptions.o: ../include/app/clOptions.hpp
app/ini.o: ../include/app/ini.hpp
app/optionparser.o: ../include/app/optionparser/optionparser.h
improc/HCIobservation.o: ../include/improc/HCIobservation.hpp ../include/improc/imageFilters.hpp
improc/ADIobservation.o: ../include/improc/ADIobservation.hpp  improc/HCIobservation.o ../include/math/gslInterpolation.hpp
improc/KLIPreduction.o: ../include/improc/KLIPreduction.hpp improc/ADIobservation.o
ioutils/fileUtils.o: ../include/ioutils/fileUtils.hpp
ioutils/fits/fitsUtils.o: ../include/ioutils/fits/fitsUtils.hpp
ioutils/fits/fitsHeaderCard.o: ../include/ioutils/fits/fitsHeaderCard.hpp ../include/ioutils/fits/fitsUtils.hpp
ioutils/fits/fitsHeader.o: ../include/ioutils/fits/fitsHeader.hpp ioutils/fits/fitsHeaderCard.o 
ipc/processInterface.o: ../include/ipc/processInterface.hpp
math/cuda/templateCublas.o: ../include/math/cuda/templateCublas.hpp
math/cuda/templateCufft.o: ../include/math/cuda/templateCufft.hpp
math/cuda/templateCurand.o: ../include/math/cuda/templateCurand.hpp
math/fit/fitGaussian.o: ../include/math/fit/fitGaussian.hpp ../include/math/fit/array2FitGaussian1D.hpp
math/fft/fftwEnvironment.o: ../include/math/fft/fftwEnvironment.hpp ../include/math/fft/fftwTemplates.hpp
math/fft/fftwTemplates.o: ../include/math/fft/fftwTemplates.hpp
math/func/bessel.o: ../include/math/func/bessel.hpp
math/func/legendre.o: ../include/math/func/legendre.hpp
math/plot/gnuPlot.o: ../include/math/plot/gnuPlot.hpp ../include/math/constants.hpp ../include/mxError.hpp ../include/mxException.hpp ../include/sys/timeUtils.hpp
sigproc/zernike.o: ../include/sigproc/zernike.hpp
sys/environment.o: ../include/sys/environment.hpp
sys/gitRepo.o: ../include/sys/gitRepo.hpp
sys/timeUtils.o: ../include/sys/timeUtils.hpp

#mxlib_uncomp_version

libmxlib.$(LIB_SUFFIX): sofa $(OBJS) $(CUDA_OBJS)
	$(MAKE) -C vendor/levmar-2.6/
ifeq ($(UNAME),Darwin)
	$(LINK.o) -dynamiclib -install_name @rpath/libmxlib.dylib -rpath @loader_path $(VENDOROBJS) $(OBJS) $(CUDA_OBJS) $(SOFA_PATH)/*.o -o libmxlib.dylib $(LDFLAGS) $(LDLIBS) $(CUDA_LIBS) $(SOFA_PATH)/libsofa_c.a
else
	$(LINK.o) -shared $(VENDOROBJS) $(OBJS) $(CUDA_OBJS) $(SOFA_PATH)/*.o -o $@ $(LDFLAGS) $(LDLIBS) $(CUDA_LIBS) 
endif

sofa:
	$(MAKE) -C $(SOFA_PATH)

install: libmxlib.$(LIB_SUFFIX)
	install -d $(LIB_PATH)
	install libmxlib.$(LIB_SUFFIX) $(LIB_PATH)

.PHONY: clean
clean:
	rm -f $(OBJS) $(CUDA_OBJS)
	$(MAKE) -C vendor/levmar-2.6/ clean
	$(MAKE) -C $(SOFA_PATH) clean

.PHONY: realclean
realclean: clean
	rm -f libmxlib.$(LIB_SUFFIX)
	rm -f $(OBJS)
	$(MAKE) -C vendor/levmar-2.6/ cleanall

#for forcing the rebuild of mxlib.o everytime, needed for git tracking
.FORCE:

# .PHONY: printvars
# printvars:
#     @$(foreach V,$(sort $(.VARIABLES)),$(if $(filter-out environment% default automatic,$(origin $V)),$(warning $V=$($V) ($(value $V)))))
