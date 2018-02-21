# Makefile for Contour Stencil Demosaicing
# Pascal Getreuer
# May 4, 2012

##
# The following three statements determine the build configuration.
# For handling different image formats, the program can be linked with
# the libjpeg, libpng, and libtiff libraries.  For each library, set
# the flags needed for linking.  To disable use of a library, comment
# its statement.  You can disable all three (BMP is always supported).
LDLIBIPOL=-lipoliio

# Radius of the mosaicked contour stencil neighborhood
RADIUS=2.5

##
# Standard make settings
SHELL=/bin/sh
CFLAGS=-O3 -ansi -pedantic -Wall -Wextra
LDFLAGS=
LDLIBS=-lm $(LDLIBIPOL)

DMCSWL1_SOURCES=dmcswl1cli.c dmcswl1.c mstencils.c \
displaycontours.c dmbilinear.c conv.c psio.c
DMBILINEAR_SOURCES=dmbilinearcli.c dmbilinear.c
MOSAIC_SOURCES=mosaic.c
IMDIFF_SOURCES=imdiff.c conv.c
GEN_MSTENCILS_SOURCES=gen_mstencils.c edge.c temsub.c

ARCHIVENAME=dmcswl1_$(shell date -u +%Y%m%d)
SOURCES=makefile.gcc makefile.vc readme.html bsd-license.txt \
inv3x3.h psio.c psio.h edge.c edge.h \
dmcswl1cli.c dmcswl1.c dmcswl1.h displaycontours.c displaycontours.h \
mosaic.c imdiff.c mstencils.tem gen_mstencils.c mstencils.c mstencils.h \
temsub.c temsub.h dmbilinearcli.c dmbilinear.c dmbilinear.h demo demo.bat \
conv.c conv.h frog.bmp doxygen.conf

##
# These statements add compiler flags to define USE_LIBJPEG, etc.,
# depending on which libraries have been specified above.
ifneq ($(LDLIBIPOL),)
	CIPOL=-DUSE_LIBJPEG -DUSE_LIBPNG -DUSE_LIBTIFF
endif

ALLCFLAGS=$(CFLAGS) $(CIPOL)
DMCSWL1_OBJECTS=$(DMCSWL1_SOURCES:.c=.o)
DMBILINEAR_OBJECTS=$(DMBILINEAR_SOURCES:.c=.o)
MOSAIC_OBJECTS=$(MOSAIC_SOURCES:.c=.o)
IMDIFF_OBJECTS=$(IMDIFF_SOURCES:.c=.o)
GEN_MSTENCILS_OBJECTS=$(GEN_MSTENCILS_SOURCES:.c=.o)
.SUFFIXES: .c .o
.PHONY: all clean rebuild srcdoc dist dist-zip

all: imintdmcswl1 imintdmbilinear imintmosaic imintdiff

imintdmcswl1: $(DMCSWL1_OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS) -s

imintdmbilinear: $(DMBILINEAR_OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS) -s

imintmosaic: $(MOSAIC_OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS) -s

imintdiff: $(IMDIFF_OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS) -s

imintgenmstencils: $(GEN_MSTENCILS_OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS) -s

mstencils.c: imintgenmstencils mstencils.tem makefile.gcc
	./imintgenmstencils $(RADIUS) mstencils.tem mstencils.c

.c.o:
	$(CC) -c $(ALLCFLAGS) $< -o $@

clean:
	$(RM) $(DMCSWL1_OBJECTS) $(DMBILINEAR_OBJECTS) \
	$(MOSAIC_OBJECTS) $(IMDIFF_OBJECTS) $(GEN_MSTENCILS_OBJECTS) \
	imintdmcswl1 imintdmbilinear imintmosaic imintdiff imintgenmstencils mstencils.c

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

