#
# Unix/Linux GCC Makefile for Levenberg - Marquardt minimization
# Under windows, use Makefile.vc for MSVC
#

CC=gcc
CONFIGFLAGS=-ULINSOLVERS_RETAIN_MEMORY


CFLAGS=$(CONFIGFLAGS) $(ARCHFLAGS) -O3 -ffast-math -funroll-loops -Wall -fPIC
LIBOBJS=lm.o Axb.o misc.o lmlec.o lmbc.o lmblec.o lmbleic.o
LIBSRCS=lm.c Axb.c misc.c lmlec.c lmbc.c lmblec.c lmbleic.c

AR=ar -r
RANLIB=ar -s

all: liblevmar.a

liblevmar.a: $(LIBOBJS)
	$(AR) liblevmar.a $(LIBOBJS)
	$(RANLIB) liblevmar.a


lm.o: lm.c lm_core.c levmar.h misc.h compiler.h
Axb.o: Axb.c Axb_core.c levmar.h misc.h
misc.o: misc.c misc_core.c levmar.h misc.h
lmlec.o: lmlec.c lmlec_core.c levmar.h misc.h
lmbc.o: lmbc.c lmbc_core.c levmar.h misc.h compiler.h
lmblec.o: lmblec.c lmblec_core.c levmar.h misc.h
lmbleic.o: lmbleic.c lmbleic_core.c levmar.h misc.h


clean:
	@rm -f $(LIBOBJS) $(DEMOBJS)

cleanall: clean
	@rm -f liblevmar.a


