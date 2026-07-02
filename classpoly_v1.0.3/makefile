CC = gcc
CFLAGS = -O3 -fomit-frame-pointer -funroll-loops -march=native -m64 -pedantic -std=gnu11
LDFLAGS ?= -static
# ff_poly is installed into a project-local prefix by the top-level build (../local).
# Override PREFIX if ff_poly lives elsewhere.
PREFIX ?= $(abspath $(CURDIR)/../local)
INCLUDES = -I$(PREFIX)/include
LIBS = -L$(PREFIX)/lib -lff_poly -lgmp -lm
INSTALL_DIR = $(PREFIX)/bin

HEADERS = bipoly.h bitmap.h class_inv.h classpoly.h classpoly_crt.h classpoly_inv.h phi_gcd.h phi_eval.h phi_fj_strings.h phi_poly.h polycosts.h velu.h \
		   bigX1.h crt.h cstd.h ecurve.h evec.h findcurve.h mpzutil.h ntutil.h prime.h iqclass.h pickprimes.h qform.h table.h tecurve.h tecurvecosts.h torcosts.h \
		   zp_poly.h zp_poly_main.h zp_poly_mul.h zp_poly_util.h zp_poly_pair.h zp2.h zp_roots.h
# zp_poly (Harvey-Sutherland large-p F_p[x]: sub-quadratic gcd/xgcd, division,
# root-finding) is folded into this library; class_inv_mpz.o (mpz_j_from_inv,
# the big-p invariant->j conversion) depends on it.
OBJECTS = bipoly.o class_inv.o classpoly.o classpoly_crt.o classpoly_inv.o compute_classpoly.o phi_eval.o phi_gcd.o phi_poly.o polycosts.o velu.o \
		  crt.o cstd.o ecurve.o ecurve_ladic.o evec.o findcurve.o mpzutil.o ntutil.o prime.o pickprimes.o iqclass.o qform.o table.o tecurve.o \
		  zp_poly_small.o zp_poly_util.o zp_poly_mul.o zp_poly_div.o zp_poly_gcd.o zp_poly_exp.o zp_poly_comp.o zp2.o zp_roots.o class_inv_mpz.o
PROGRAMS = classpoly invtoj

all: $(PROGRAMS)

install: $(PROGRAMS)
	mkdir -p $(INSTALL_DIR)
	cp $(PROGRAMS) $(INSTALL_DIR)

clean:
	rm -f *.o
	rm -f $(PROGRAMS)

##### executables

classpoly: classpoly_program.o $(OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $< $(OBJECTS) $(LIBS)

invtoj: invtoj.o $(OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $< $(OBJECTS) $(LIBS)

##### C modules (one generic rule covers every object file)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<
