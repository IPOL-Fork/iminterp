# GCC makefile for sinterp, Contour stencil interpolation demo
# Pascal Getreuer
# May 4, 2012

##
# The following three statements determine the build configuration.
# For handling different image formats, the program can be linked with
# the libjpeg, libpng, and libtiff libraries.  For each library, set
# the flags needed for linking.  To disable use of a library, comment
# its statement.  You can disable all three (BMP is always supported).
LDLIBIPOL=-lipoliio

##
# Standard make settings
SHELL=/bin/sh
CFLAGS=-O3 -ansi -pedantic -Wall -Wextra
LDFLAGS=
LDLIBS=-lm $(LDLIBIPOL)

##
# These statements add compiler flags to define USE_LIBJPEG, etc.,
# depending on which libraries have been specified above.
ifneq ($(LDLIBIPOL),)
	CIPOL=-DUSE_LIBIPOL -DUSE_LIBJPEG -DUSE_LIBPNG -DUSE_LIBTIFF
endif

ALLCFLAGS=$(CFLAGS) $(CIPOL)

SINTERP_SOURCES=sinterpcli.c sinterp.c sset.c invmat.c svd2x2.c pen.c
IMCOARSEN_SOURCES=imcoarsen.c
IMDIFF_SOURCES=imdiff.c conv.c
NNINTERP_SOURCES=nninterpcli.c nninterp.c

ARCHIVENAME=sinterp_$(shell date -u +%Y%m%d)
SOURCES=conv.c conv.h imcoarsen.c imdiff.c \
invmat.c invmat.h nninterp.c nninterp.h nninterpcli.c pen.c pen.h svd2x2.c svd2x2.h sinterp.c \
sinterp.h sinterpcli.c sset.c sset.h readme.html bsd-license.txt \
makefile.gcc makefile.vc doxygen.conf demo demo.bat frog-hr.bmp
SINTERP_OBJECTS=$(SINTERP_SOURCES:.c=.o)
IMCOARSEN_OBJECTS=$(IMCOARSEN_SOURCES:.c=.o)
IMDIFF_OBJECTS=$(IMDIFF_SOURCES:.c=.o)
NNINTERP_OBJECTS=$(NNINTERP_SOURCES:.c=.o)
.SUFFIXES: .c .o
.PHONY: all clean rebuild srcdoc dist dist-zip

all: iminterps imintcoarsen imintdiff iminterpnn

iminterps: $(SINTERP_OBJECTS)
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
	$(RM) $(SINTERP_OBJECTS) $(IMCOARSEN_OBJECTS) \
	$(IMDIFF_OBJECTS) $(NNINTERP_OBJECTS) \
	iminterps imintcoarsen imintdiff iminterpnn

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
