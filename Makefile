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

-include mk/Common.mk

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
                 sys \
                 wfp \
                 ioutils \
                 mxException.hpp \
                 mxError.hpp \
                 mxlib.hpp\
                 mxlib_uncomp_version.h
                
all: lib

.PHONY: mxlib_uncomp_version
mxlib_uncomp_version:
	@bash ./gengithead.sh ./ ./include/mxlib_uncomp_version.h MXLIB_UNCOMP

.PHONY: mxlib_comp_version
mxlib_comp_version:
	@bash ./gengithead.sh ./ ./source/mxlib_comp_version.h MXLIB_COMP
	
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

lib: mxlib_uncomp_version mxlib_comp_version
	cd source; ${MAKE}
	
install: all mxlib_uncomp_version
	cd source; ${MAKE} install
	install -d $(INCLUDE_PATH)/mx
	install -d $(BIN_PATH)
	install gengithead.sh $(BIN_PATH)/
	cp -r include/* $(INCLUDE_PATH)/mx/

.PHONY: clean
clean:
	rm -f *.o *~
	rm -f include/mxlib_uncomp_version.h
	rm -f include/mxlib_comp_version.h
	$(MAKE) -C source clean

.PHONY: realclean
realclean: clean
	$(MAKE) -C source realclean
