INCLUDE_PATH = $(HOME)/include
LIB_PATH = $(HOME)/lib
MXLIB_EXLIBS = -lsofa_c -L/usr/lib64/ -lrt -lcfitsio -L/usr/local/lib -lfftw3 -lfftw3f -lfftw3_threads -lfftw3f_threads -lboost_system -lboost_filesystem -L/usr/local/atlas/lib/ -ltatlas  -lgsl -lgslcblas 

OPTIMIZE = -O3 -fopenmp 
#-ffast-math -Wno-ignored-attributes

CXX = g++

CFLAGS += --std=c99 -D_XOPEN_SOURCE=600  -fPIC 
#-fno-diagnostics-show-caret
CXXFLAGS += --std=c++0x -D_XOPEN_SOURCE=600 -fPIC 
#-fno-diagnostics-show-caret

INCLUDE = -I../include -I$(INCLUDE_PATH) -I/usr/local/atlas/include 

.c.o:
	$(CC) $(OPTIMIZE) $(CFLAGS) $(INCLUDE) -c $<

.cpp.o:
	$(CXX) $(OPTIMIZE) $(CXXFLAGS) $(INCLUDE) -c $<

# programs to be made
TARGET = $(targname)

OBJS = $(targname).o


all: git_version $(OBJS) 
	$(CXX) -o ../bin/$(TARGET)  $(OBJS) $(OPTIMIZE) $(CXXFLAGS) -L$(LIB_PATH) -lmxlib  $(MXLIB_EXLIBS)

pch: 
	$(CXX) $(OPTIMIZE) $(CXXFLAGS) $(INCLUDE) -x c++-header ../include/mxAOAnalytic.hpp

.PHONY: git_version
git_version:
	gengithead.sh ../ ../include/mxaoa_version.h MXAOANALYTIC


.PHONEY: clean
clean:
	rm $(TARGET)
	rm -f *.o 
	rm -f *~

