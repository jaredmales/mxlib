##########################################
##                                      ##
##          Makefile for mxlib          ##
##                                      ##
##########################################

# Are you here to customize the build?
#
# Option 1: Override the options on the command line
#     e.g. `make install PREFIX=$HOME`
# Option 2: `make setup`
#     This creates ./local/{Common,MxLib,MxApp}.mk and lists options.
#     Define those options in the ./local Makefiles. These files are
#     *not* tracked in git, so your tweaks won't cause warnings about
#     uncommitted changes.

include mk/Common.mk
include mk/MxLib.mk


INC_TO_INSTALL = ao \
                 app \
                 astro \
                 cuda \
                 fft \
                 improc \
                 ipc \
                 math \
                 meta \
                 sigproc \
                 wfp \
                 ioutils \
                 eigenUtils.hpp \
                 environment.hpp \
                 gnuPlot.hpp \
                 gslInterpolation.hpp \
                 imagingArray.hpp \
                 mxException.hpp \
                 mxError.hpp \
                 mxlib.hpp\
                 mxlib_uncomp_version.h\
                 ompLoopWatcher.hpp \
                 timeUtils.hpp 

all: install

.PHONY: mxlib_uncomp_version
mxlib_uncomp_version:
	@bash ./gengithead.sh ./ ./include/mxlib_uncomp_version.h MXLIB_UNCOMP

.PHONY: setup
setup:
	@for file in ./local/*.example.mk; do \
		dest=$$(echo $$file | sed 's/.example//'); \
		if [ ! -e $$dest ]; then cp -v $$file $$dest; fi \
	done
	@echo "***\nBuild settings available in local/Common.mk\n***"
	@grep "?=" mk/Common.mk || true
	@echo "***"
	@echo "Build settings available in local/MxLib.mk\n***"
	@grep "?=" mk/MxLib.mk || true
	@echo "***"
	@echo "Build settings available in local/MxApp.mk\n***"
	@grep  "?=" mk/MxApp.mk || true
	@echo "***"


install: mxlib_uncomp_version
	install -d $(INCLUDE_PATH)/mx
	install gengithead.sh $(BIN_PATH)
	for file in ${INC_TO_INSTALL}; do \
	  (cp -r include/$$file $(INCLUDE_PATH)/mx) || break; \
	done

clean:
	rm -f *.o *~
	rm -f include/mxlib_uncomp_version.h
