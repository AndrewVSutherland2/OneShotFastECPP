CC = gcc
# class_inv.h (invariant codes + inv_string) and ff_poly.h / zp_poly.h are needed
# by cminv/cm_method; both libraries are staged in ../local by the top-level make.
CPDIR = ../classpoly_v1.0.3
CPINC = -isystem $(CPDIR) -isystem ../local/include
CFLAGS = -O3 -march=native -m64 -funroll-loops -std=gnu11 -Wall -Wextra -fopenmp $(CPINC)
LIBS = -lgmp -lm

# cminv's production path calls classpoly's mpz_j_from_inv (class_inv_mpz.c), which
# pulls in classpoly's phi-file machinery and the zp_poly code (now part of the
# classpoly library).  Link everything classpoly has except its two program mains.
CPOBJS = $(filter-out $(CPDIR)/classpoly_program.o $(CPDIR)/invtoj.o, $(wildcard $(CPDIR)/*.o))
CPLIBS = -L../local/lib -lff_poly $(LIBS)
# the zp_poly objects fproot's sub-quadratic gcd path needs (subset of CPOBJS)
ZPOBJS = $(addprefix $(CPDIR)/, zp_poly_small.o zp_poly_util.o zp_poly_mul.o zp_poly_div.o zp_poly_gcd.o zp_poly_exp.o zp_poly_comp.o zp2.o zp_roots.o)

PROGRAMS = cmsearch dscan smoothtest roottest hdroot cm_method oneshot
TESTPROGS = invjtest cmjtest asmtest zpbench
GENERIC  = cmsearch dscan smoothtest roottest hdroot oneshot invjtest asmtest

all: $(PROGRAMS)
tests: $(TESTPROGS)

cmsearch:   cmsearch.o cornacchia.o
dscan:      dscan.o cornacchia.o
smoothtest: smoothtest.o smooth.o
roottest:   roottest.o fproot.o
hdroot:     hdroot.o fproot.o
oneshot:    oneshot.o cornacchia.o smooth.o fproot.o curve.o
invjtest:   invjtest.o invj.o fproot.o
asmtest:    asmtest.o curve.o fproot.o cornacchia.o smooth.o

$(GENERIC):
	$(CC) $(CFLAGS) -o $@ $^ $(ZPOBJS) $(LIBS)

cm_method:  cm_method.o cornacchia.o fproot.o cminv.o invj.o
	$(CC) $(CFLAGS) -o $@ $^ $(CPOBJS) $(CPLIBS)

cmjtest:    cmjtest.o cminv.o invj.o fproot.o
	$(CC) $(CFLAGS) -o $@ $^ $(CPOBJS) $(CPLIBS)

zpbench:    zpbench.o fproot.o
	$(CC) $(CFLAGS) -o $@ $^ $(ZPOBJS) $(LIBS)

# header dependencies
cornacchia.o cmsearch.o dscan.o: cornacchia.h
smooth.o smoothtest.o:           smooth.h
fproot.o roottest.o hdroot.o:    fproot.h
invj.o cminv.o invjtest.o cmjtest.o: invj.h fproot.h
curve.o asmtest.o:               curve.h fproot.h cornacchia.h smooth.h
cminv.o cm_method.o:             cminv.h
oneshot.o:                       cornacchia.h smooth.h fproot.h curve.h

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o $(PROGRAMS) $(TESTPROGS)
