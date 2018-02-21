# GCC makefile for TV denoising IPOL demo
# Pascal Getreuer
# May 16, 2012

# The following three statements determine the build configuration.
# For handling different image formats, the program can be linked with
# the libjpeg, libpng, and libtiff libraries.  For each library, set
# the flags needed for linking.  To disable use of a library, comment
# its statement.  You can disable all three (BMP is always supported).
LDLIBIPOL=-lipoliio

TVREG_FLAGS=-DTVREG_DENOISE -DTVREG_NONGAUSSIAN -DNUM_SINGLE

##
# Standard make settings
CFLAGS=-O3 -ansi -pedantic -Wall -Wextra $(TVREG_FLAGS)
LDFLAGS=
LDLIB=-lm $(LDLIBIPOL)

TVDENOISE_SOURCES=tvdenoise.c tvreg.c
IMNOISE_SOURCES=imnoise.c randmt.c
IMDIFF_SOURCES=imdiff.c conv.c

ARCHIVENAME=tvdenoise_$(shell date -u +%Y%m%d)
SOURCES=tvdenoise.c tvreg.c tvreg.h tvregopt.h dsolve_inc.c zsolve_inc.c \
usolve_gs_inc.c imnoise.c randmt.c randmt.h imdiff.c conv.c conv.h \
num.h makefile.gcc makefile.vc \
readme.txt code_overview.txt license.txt doxygen.conf einstein.bmp example.sh

## 
# These statements add compiler flags to define USE_LIBJPEG, etc.,
# depending on which libraries have been specified above.
ifneq ($(LDLIBIPOL),)
	CIPOL=-DUSE_LIBJPEG -DUSE_LIBPNG -DUSE_LIBTIFF
endif

ALLCFLAGS=$(CFLAGS) $(CIPOL)
TVDENOISE_OBJECTS=$(TVDENOISE_SOURCES:.c=.o)
IMNOISE_OBJECTS=$(IMNOISE_SOURCES:.c=.o)
IMDIFF_OBJECTS=$(IMDIFF_SOURCES:.c=.o)
.SUFFIXES: .c .o
.PHONY: all clean rebuild srcdoc dist dist-zip

all: iminttvdenoise imintnoise imintdiff

iminttvdenoise: $(TVDENOISE_OBJECTS)
	$(CC) $(LDFLAGS) $^ $(LDLIB) -o $@ -s

imintnoise: $(IMNOISE_OBJECTS)
	$(CC) $(LDFLAGS) $^ $(LDLIB) -o $@ -s

imintdiff: $(IMDIFF_OBJECTS)
	$(CC) $(LDFLAGS) $^ $(LDLIB) -o $@ -s

.c.o:
	$(CC) -c $(ALLCFLAGS) $< -o $@

clean:
	$(RM) $(TVDENOISE_OBJECTS) $(IMNOISE_OBJECTS) \
	$(IMDIFF_OBJECTS) iminttvdenoise imintnoise imintdiff

rebuild: clean all

srcdoc:
	doxygen doxygen.conf

dist: $(SOURCES)
	-rm -rf $(ARCHIVENAME)
	mkdir $(ARCHIVENAME)
	ln $(SOURCES) $(ARCHIVENAME)
	tar vchzf $(ARCHIVENAME).tar.gz $(ARCHIVENAME)
	-rm -rf $(ARCHIVENAME)

dist-zip: $(SOURCES)
	-rm -rf $(ARCHIVENAME)
	mkdir $(ARCHIVENAME)
	ln $(SOURCES) $(ARCHIVENAME)
	-rm -f $(ARCHIVENAME).zip
	zip -r9 $(ARCHIVENAME).zip $(ARCHIVENAME)/*
	-rm -rf $(ARCHIVENAME)
