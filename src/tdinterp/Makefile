# GCC makefile for tdinterp, Tensor-driven diffusion interpolation demo
# Pascal Getreuer
# May 4, 2012

# The FFTW3 single precision library (http://www.fftw.org) is required.
# Set the flags needed for linking.
LDFFTW3=-lfftw3f

# The following three statements determine the build configuration.
# For handling different image formats, the program can be linked with
# the libjpeg, libpng, and libtiff libraries.  For each library, set
# the flags needed for linking.  To disable use of a library, comment
# its statement.  You can disable all three (BMP is always supported).
LDLIBIPOL=-lipoliio

##
# Set this line to compile with OpenMP multithreading.  Comment the
# line to disable OpenMP.
OPENMP=-fopenmp

##
# Standard make settings
SHELL=/bin/sh
CFLAGS=-O3 -ansi -pedantic -Wall -Wextra $(OPENMP)
LDFLAGS=
LDLIBS=-lm $(LDFFTW3) $(LDLIBIPOL) $(OPENMP)

##
# These statements add compiler flags to define USE_LIBJPEG, etc.,
# depending on which libraries have been specified above.
ifneq ($(LDLIBIPOL),)
	CIPOL=-DUSE_LIBJPEG -DUSE_LIBPNG -DUSE_LIBTIFF
endif

ALLCFLAGS=$(CFLAGS) $(CIPOL)

TDINTERP_SOURCES=tdinterpcli.c tdinterp.c conv.c finterp.c
IMCOARSEN_SOURCES=imcoarsen.c
IMDIFF_SOURCES=imdiff.c conv.c
NNINTERP_SOURCES=nninterpcli.c nninterp.c

ARCHIVENAME=tdinterp_$(shell date -u +%Y%m%d)
SOURCES=conv.c conv.h imcoarsen.c imdiff.c \
tdinterpcli.c tdinterp.c tdinterp.h finterp.c finterp.h \
nninterpcli.c nninterp.c nninterp.h demo demo.bat frog-hr.bmp \
readme.html bsd-license.txt makefile.gcc makefile.vc doxygen.conf
TDINTERP_OBJECTS=$(TDINTERP_SOURCES:.c=.o)
IMCOARSEN_OBJECTS=$(IMCOARSEN_SOURCES:.c=.o)
IMDIFF_OBJECTS=$(IMDIFF_SOURCES:.c=.o)
NNINTERP_OBJECTS=$(NNINTERP_SOURCES:.c=.o)
.SUFFIXES: .c .o
.PHONY: all clean rebuild srcdoc dist dist-zip

all: iminterptd imintcoarsen imintdiff iminterpnn

iminterptd: $(TDINTERP_OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS) -s

imintcoarsen: $(IMCOARSEN_OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS) -s

imintdiff: $(IMDIFF_OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS) -s

iminterpnn: $(NNINTERP_OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS) -s

.c.o:
	$(CC) -c $(ALLCFLAGS) $< -o $@

clean:
	$(RM) $(TDINTERP_OBJECTS) $(IMCOARSEN_OBJECTS) $(IMDIFF_OBJECTS) $(NNINTERP_OBJECTS) \
	iminterptd imintcoarsen imintdiff iminterpnn

rebuild: clean all

# Source documentation with Doxygen
srcdoc: $(SOURCES)
	doxygen doxygen.conf

dist: $(SOURCES)
	-rm -rf $(ARCHIVENAME)
	mkdir $(ARCHIVENAME)
	ln $(SOURCES) $(ARCHIVENAME)
	tar vchzf $(ARCHIVENAME).tgz $(ARCHIVENAME)
	-rm -rf $(ARCHIVENAME)

dist-zip: $(SOURCES)
	-rm -rf $(ARCHIVENAME)
	mkdir $(ARCHIVENAME)
	ln $(SOURCES) $(ARCHIVENAME)
	-rm -f $(ARCHIVENAME).zip
	zip -r9 $(ARCHIVENAME).zip $(ARCHIVENAME)/*
	-rm -rf $(ARCHIVENAME)
