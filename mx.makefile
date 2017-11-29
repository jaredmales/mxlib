########################################################################
## This is makefile for building mxlib-based single-file applications ##
########################################################################

# Instructions:
# 1) define the environment variable MXMAKEFILEINC to point to mx.makefile.inc 
# 2) define the environment variable MXMAKEFILE to point to this file.
# 3) invoke it make with make -B -f $MXMAKEFILE t=source
#        -- here source refers to source.cpp, which is compiled to an executable named source 
# 4) run it with ./source


include $(MXMAKEFILEINC)

# programs to be made
TARGET = $(t)

OBJS = $(t).o


all: $(OBJS) 
	$(CXX) -o $(TARGET)  $(OBJS) $(OPTIMIZE) $(CXXFLAGS) -L$(LIB_PATH) -lmxlib  $(MXLIB_EXLIBS)

install: all
	install -d $(BIN_PATH)
	install $(TARGET) $(BIN_PATH)
	
.PHONEY: clean
clean:
	rm $(TARGET)
	rm -f *.o 
	rm -f *~

