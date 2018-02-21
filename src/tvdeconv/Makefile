# GCC makefile for TV deconvolution IPOL demo
# Pascal Getreuer
# May 16, 2012

# The FFTW3 library (http://www.fftw.org) is required for deconvolution.
# Set the flags needed for linking.
LDLIBFFTW3=-lfftw3 -lfftw3f

# The following three statements determine the build configuration.
# For handling different image formats, the program can be linked with
# the libjpeg, libpng, and libtiff libraries.  For each library, set
# the flags needed for linking.  To disable use of a library, comment
# its statement.  You can disable all three (BMP is always supported).
LDLIBIPOL=-lipoliio

TVREG_FLAGS=-DTVREG_DECONV -DTVREG_NONGAUSSIAN -DNUM_SINGLE

##
# Standard make settings
CFLAGS=-O3 -ansi -pedantic -Wall -Wextra $(TVREG_FLAGS)
LDFLAGS=
LDLIB=-lm $(LDLIBFFTW3) $(LDLIBIPOL)

#CC=gcc

TVDECONV_SOURCES=tvdeconv.c tvreg.c kernels.c cliio.c
IMBLUR_SOURCES=imblur.c kernels.c randmt.c cliio.c
IMDIFF_SOURCES=imdiff.c conv.c

ARCHIVENAME=tvdeconv_$(shell date -u +%Y%m%d)
SOURCES=tvdeconv.c tvreg.c tvreg.h tvregopt.h dsolve_inc.c zsolve_inc.c \
usolve_dct_inc.c usolve_dft_inc.c util_deconv.h imblur.c randmt.c randmt.h \
imdiff.c conv.c conv.h kernels.c kernels.h cliio.c cliio.h \
num.h makefile.gcc makefile.vc \
readme.txt license.txt doxygen.conf einstein.bmp example.sh

##
# These statements add compiler flags to define USE_LIBJPEG, etc.,
# depending on which libraries have been specified above.
ifneq ($(LDLIBIPOL),)
	CIPOL=-DUSE_LIBJPEG -DUSE_LIBPNG -DUSE_LIBTIFF
endif

ALLCFLAGS=$(CFLAGS) $(CIPOL)
TVDECONV_OBJECTS=$(TVDECONV_SOURCES:.c=.o)
IMBLUR_OBJECTS=$(IMBLUR_SOURCES:.c=.o)
IMDIFF_OBJECTS=$(IMDIFF_SOURCES:.c=.o)
.SUFFIXES: .c .o
.PHONY: all clean rebuild srcdoc dist dist-zip

all: iminttvdeconv imintblur imintdiff

iminttvdeconv: $(TVDECONV_OBJECTS)
	$(CC) $(LDFLAGS) $^ $(LDLIB) -o $@ -s

imintblur: $(IMBLUR_OBJECTS)
	$(CC) $(LDFLAGS) $^ $(LDLIB) -o $@ -s

imintdiff: $(IMDIFF_OBJECTS)
	$(CC) $(LDFLAGS) $^ $(LDLIB) -o $@ -s

.c.o:
	$(CC) -c $(ALLCFLAGS) $< -o $@

clean:
	$(RM) $(TVDECONV_OBJECTS) $(IMBLUR_OBJECTS) \
	$(IMDIFF_OBJECTS) iminttvdeconv imintblur imintdiff

rebuild: clean all

srcdoc:
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
