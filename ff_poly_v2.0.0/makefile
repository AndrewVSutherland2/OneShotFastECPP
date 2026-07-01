CC = gcc
CFLAGS = -O3 -fomit-frame-pointer -funroll-loops -march=native -m64 -pedantic -std=gnu11
LDFLAGS = -static
INCLUDES = -I.
LIBS = -lgmp -lm
# Install into a project-local prefix by default so nothing touches /usr/local.
# Override with `make PREFIX=/some/where install` if desired.
PREFIX ?= $(abspath $(CURDIR)/../local)
INSTALL_ROOT = $(PREFIX)

HEADERS = asm.h ff.h fff.h ffmontgomery64.h ff2k.h ffext.h ffpoly.h ffpolyalloc.h ffpolybig.h ffpolyfromroots.h ffpolymultipoint.h ffpolysmall.h cstd.h ntutil.h polyparse.h
OBJECTS = ff.o ff2k.o  ffext.o ffpoly.o ffpolyalloc.o ffpolybig.o ffpolyfromroots.o ffpolymultipoint.o ffpolysmall.o polyparse.o

ZNPHEADERS = zn_poly/wide_arith.h zn_poly/zn_poly.h zn_poly/zn_poly_internal.h
ZNPOBJECTS = zn_poly/array.o zn_poly/ks_support.o zn_poly/misc.o zn_poly/mpn_mulmid.o zn_poly/mul.o zn_poly/mul_fft.o zn_poly/mul_ks.o zn_poly/mulmid.o zn_poly/mulmid_ks.o zn_poly/nuss.o zn_poly/pack.o zn_poly/pmf.o zn_poly/pmfvec_fft.o zn_poly/tuning.o zn_poly/zn_mod.o

all: libff_poly.a

clean:
	rm -f *.o
	rm -f zn_poly/*.o
	rm -f libff_poly.a
	
install: all
	mkdir -p $(INSTALL_ROOT)/include/ff_poly
	mkdir -p $(INSTALL_ROOT)/lib
	cp -v ff_poly.h $(INSTALL_ROOT)/include
	cp -v $(HEADERS) $(INSTALL_ROOT)/include/ff_poly
	cp -v libff_poly.a $(INSTALL_ROOT)/lib
	
##### ff_poly library

libff_poly.a: $(OBJECTS) $(ZNPOBJECTS)
	ar -r libff_poly.a $(OBJECTS) $(ZNPOBJECTS)
	ranlib libff_poly.a
	
##### executables

##### ff_poly modules

ff.o: ff.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

ff2k.o : ff2k.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

ffext.o: ffext.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

ffpolyalloc.o: ffpolyalloc.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

ffpoly.o: ffpoly.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

ffpolybig.o: ffpolybig.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

ffpolyfromroots.o: ffpolyfromroots.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

ffpolymultipoint.o: ffpolymultipoint.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

ffpolysmall.o: ffpolysmall.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

polyparse.o: polyparse.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

ff_poly.h: $(HEADERS)
	touch ff_poly.h

##### zn_poly modules

zn_poly/array.o: zn_poly/array.c $(ZNPOLYHEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

zn_poly/ks_support.o: zn_poly/ks_support.c $(ZNPOLYHEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

zn_poly/misc.o: zn_poly/misc.c $(ZNPOLYHEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

zn_poly/mpn_mulmid.o: zn_poly/mpn_mulmid.c $(ZNPOLYHEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

zn_poly/mul.o: zn_poly/mul.c $(ZNPOLYHEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

zn_poly/mul_fft.o: zn_poly/mul_fft.c $(ZNPOLYHEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

zn_poly/mul_ks.o: zn_poly/mul_ks.c $(ZNPOLYHEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

zn_poly/mulmid.o: zn_poly/mulmid.c $(ZNPOLYHEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

zn_poly/mulmid_ks.o: zn_poly/mulmid_ks.c $(ZNPOLYHEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

zn_poly/nuss.o: zn_poly/nuss.c $(ZNPOLYHEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

zn_poly/pack.o: zn_poly/pack.c $(ZNPOLYHEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

zn_poly/pmf.o: zn_poly/pmf.c $(ZNPOLYHEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

zn_poly/pmfvec_fft.o: zn_poly/pmfvec_fft.c $(ZNPOLYHEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

zn_poly/zn_mod.o: zn_poly/zn_mod.c $(ZNPOLYHEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

zn_poly/tuning.o: zn_poly/tuning.c $(ZNPOLYHEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<
