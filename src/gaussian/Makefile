# GCC makefile for Gaussian Convolution
# Pascal Getreuer
# April 28, 2013

# The FFTW3 library (http://www.fftw.org) is required.
# Set the flags needed for linking with fftw3 and fftw3f.
LDFFTW3=-lfftw3 -lfftw3f

# The following three statements determine the build configuration.
# For handling different image formats, the program can be linked with
# the libjpeg, libpng, and libtiff libraries.  For each library, set
# the flags needed for linking.  To disable use of a library, comment
# its statement.  You can disable all three (BMP is always supported).
LDLIBIPOL=-lipoliio

# Set this line to compute in single precision instead of double precision
# NUM=-DNUM_SINGLE

# Make settings
SHELL=/bin/sh
CFLAGS=-O3 -ansi -pedantic -Wall -Wextra $(NUM)
LDFLAGS=$(OPENMP)
LDLIBS=-lm $(LDFFTW3) $(LDLIBIPOL)

GAUSSIAN_BENCH_SOURCES=gaussian_bench.c strategy_gaussian_conv.c \
gaussian_conv_fir.c gaussian_conv_dct.c gaussian_conv_am.c \
gaussian_conv_deriche.c gaussian_conv_vyv.c \
gaussian_conv_box.c gaussian_conv_ebox.c \
gaussian_conv_sii.c gaussian_short_conv.c \
filter_util.c erfc_cody.c inverfc_acklam.c invert_matrix.c
GAUSSIAN_DEMO_SOURCES=gaussian_demo.c \
gaussian_conv_fir.c gaussian_conv_dct.c \
gaussian_conv_am.c gaussian_conv_deriche.c \
gaussian_conv_vyv.c gaussian_conv_box.c gaussian_conv_ebox.c \
gaussian_conv_sii.c gaussian_short_conv.c \
filter_util.c erfc_cody.c inverfc_acklam.c invert_matrix.c
IMDIFF_SOURCES=imdiff.c

ARCHIVENAME=gaussian_$(shell date -u +%Y%m%d)
SOURCES=gaussian_bench.c gaussian_demo.c \
gaussian_conv_fir.c gaussian_conv_fir.h \
gaussian_conv_dct.c gaussian_conv_dct.h \
gaussian_conv_am.c gaussian_conv_am.h \
gaussian_conv_deriche.c gaussian_conv_deriche.h \
gaussian_conv_vyv.c gaussian_conv_vyv.h \
gaussian_conv_box.c gaussian_conv_box.h \
gaussian_conv_ebox.c gaussian_conv_ebox.h \
gaussian_conv_sii.c gaussian_conv_sii.h \
gaussian_short_conv.c gaussian_short_conv.h \
strategy_gaussian_conv.c strategy_gaussian_conv.h \
filter_util.c filter_util.h num.h complex_arith.h \
erfc_cody.c erfc_cody.h inverfc_acklam.c inverfc_acklam.h \
invert_matrix.c invert_matrix.h imdiff.c makefile.gcc \
README.txt demo.sh bench.sh plotimpulse.gp einstein.png \
doxygen.conf doxygenlayout.xml BSD_simplified.txt GPLv3.txt

# These statements add compiler flags to define USE_LIBJPEG, etc.,
# depending on which libraries have been specified above.
ifneq ($(LDLIBIPOL),)
	CIPOL=-DUSE_LIBJPEG -DUSE_LIBPNG -DUSE_LIBTIFF
endif

ALLCFLAGS=$(CFLAGS) $(CIPOL)
GAUSSIAN_BENCH_OBJECTS=$(GAUSSIAN_BENCH_SOURCES:.c=.o)
GAUSSIAN_DEMO_OBJECTS=$(GAUSSIAN_DEMO_SOURCES:.c=.o)
IMDIFF_OBJECTS=$(IMDIFF_SOURCES:.c=.o)
.SUFFIXES: .c .o
.PHONY: all clean rebuild srcdoc dist dist-zip dist-xz

all: imintgaussianbench imintgaussianconv imintdiff

imintgaussianbench: $(GAUSSIAN_BENCH_OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS) -s

imintgaussianconv: $(GAUSSIAN_DEMO_OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS) -s

imintdiff: $(IMDIFF_OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS) -s

.c.o:
	$(CC) -c $(ALLCFLAGS) $< -o $@

clean:
	$(RM) $(GAUSSIAN_BENCH_OBJECTS) $(GAUSSIAN_DEMO_OBJECTS) $(IMDIFF_OBJECTS) \
	imintgaussianbench imintgaussianconv imintdiff

rebuild: clean all

# Source documentation with Doxygen
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
	-rm $(ARCHIVENAME).zip
	zip -r9 $(ARCHIVENAME).zip $(ARCHIVENAME)/*
	-rm -rf $(ARCHIVENAME)

dist-xz: $(SOURCES)
	-rm -rf $(ARCHIVENAME)
	mkdir $(ARCHIVENAME)
	ln $(SOURCES) $(ARCHIVENAME)
	tar vchJf $(ARCHIVENAME).tar.xz $(ARCHIVENAME)
	-rm -rf $(ARCHIVENAME)
