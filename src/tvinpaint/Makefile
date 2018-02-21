# GCC makefile for TV inpainting IPOL demo
# Pascal Getreuer
# May 16, 2012

# The following three statements determine the build configuration.
# For handling different image formats, the program can be linked with
# the libjpeg, libpng, and libtiff libraries.  For each library, set
# the flags needed for linking.  To disable use of a library, comment
# its statement.  You can disable all three (BMP is always supported).
LDLIBIPOL=-lipoliio

TVREG_FLAGS=-DTVREG_INPAINT -DNUM_SINGLE

##
# Standard make settings
CFLAGS=-O3 -ansi -pedantic -Wall -Wextra $(TVREG_FLAGS)
LDFLAGS=
LDLIB=-lm $(LDLIBIPOL)

TVINPAINT_SOURCES=tvinpaint.c tvreg.c
RANDMASK_SOURCES=randmask.c randmt.c drawtext.c
APPLYMASK_SOURCES=applymask.c

ARCHIVENAME=tvinpaint_$(shell date -u +%Y%m%d)
SOURCES=tvinpaint.c tvreg.c tvreg.h tvregopt.h dsolve_inc.c usolve_gs_inc.c \
num.h randmask.c randmt.c randmt.h drawtext.c drawtext.h applymask.c \
makefile.gcc makefile.vc readme.txt \
code_overview.txt BSD_simplified.txt GPLv3.txt doxygen.conf mountain.bmp \
example.sh

## 
# These statements add compiler flags to define USE_LIBJPEG, etc.,
# depending on which libraries have been specified above.
ifneq ($(LDLIBIPOL),)
	CIPOL=-DUSE_LIBJPEG -DUSE_LIBPNG -DUSE_LIBTIFF
endif

ALLCFLAGS=$(CFLAGS) $(CIPOL)
TVINPAINT_OBJECTS=$(TVINPAINT_SOURCES:.c=.o)
RANDMASK_OBJECTS=$(RANDMASK_SOURCES:.c=.o)
APPLYMASK_OBJECTS=$(APPLYMASK_SOURCES:.c=.o)
.SUFFIXES: .c .o
.PHONY: all clean rebuild srcdoc dist dist-zip

all: iminttvinpaint imintmaskrand imintmaskapply

iminttvinpaint: $(TVINPAINT_OBJECTS)
	$(CC) $(LDFLAGS) $^ $(LDLIB) -o $@ -s

imintmaskrand: $(RANDMASK_OBJECTS)
	$(CC) $(LDFLAGS) $^ $(LDLIB) -o $@ -s

imintmaskapply: $(APPLYMASK_OBJECTS)
	$(CC) $(LDFLAGS) $^ $(LDLIB) -o $@ -s

.c.o:
	$(CC) -c $(ALLCFLAGS) $< -o $@

clean:
	$(RM) $(TVINPAINT_OBJECTS) $(RANDMASK_OBJECTS) $(APPLYMASK_OBJECTS) \
	iminttvinpaint imintmaskrand imintmaskapply

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
