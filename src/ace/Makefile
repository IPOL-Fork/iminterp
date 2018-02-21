# Makefile for ACE
# Pascal Getreuer
# October 28, 2012

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
# Make settings
SHELL=/bin/sh
CFLAGS=-O3 -ansi -pedantic -Wall -Wextra $(OPENMP)
LDFLAGS=$(OPENMP)
LDLIBS=-lm $(LDFFTW3) $(LDLIBIPOL)

##
# These statements add compiler flags to define USE_LIBJPEG, etc.,
# depending on which libraries have been specified above.
ifneq ($(LDLIBIPOL),)
	CIPOL=-DUSE_LIBJPEG -DUSE_LIBPNG -DUSE_LIBTIFF
endif

ALLCFLAGS=$(CFLAGS) $(CIPOL)

ACE_SOURCES=acecli.c ace.c
HISTEQ_SOURCES=histeq.c

ARCHIVENAME=ace_$(shell date -u +%Y%m%d)
SOURCES=acecli.c ace.c ace.h slopecoeff_inc.c histeq.c \
makefile.gcc doxygen.conf \
BSD_simplified.txt GPLv3.txt avs.jpg example.sh readme.txt
ACE_OBJECTS=$(ACE_SOURCES:.c=.o)
HISTEQ_OBJECTS=$(HISTEQ_SOURCES:.c=.o)

.SUFFIXES: .c .o
.PHONY: all clean rebuild srcdoc dist dist-zip

all: imintace iminthisteq

imintace: $(ACE_OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS) -s

iminthisteq: $(HISTEQ_OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS) -s

.c.o:
	$(CC) -c $(ALLCFLAGS) $< -o $@

clean:
	$(RM) $(ACE_OBJECTS) $(HISTEQ_OBJECTS) imintace iminthisteq

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
