CC = gcc
# class_inv.h (invariant codes + inv_string) and ff_poly.h / zp_poly.h are needed
# by cminv/cm_method; both libraries are staged in ../local by the top-level make.
CPDIR = ../classpoly_v1.0.3
CPINC = -isystem $(CPDIR) -isystem ../local/include
CFLAGS = -O3 -march=native -m64 -funroll-loops -std=gnu11 -Wall -Wextra -fopenmp $(CPINC)
LIBS = -L../local/lib -lzp_poly -lgmp -lm

# cminv's production path calls classpoly's mpz_j_from_inv (class_inv_mpz.c), which
# pulls in classpoly's phi-file machinery and the zp_poly library.  Link everything
# classpoly has except its two program mains (mirrors the classpoly binary link).
CPOBJS = $(filter-out $(CPDIR)/classpoly_program.o $(CPDIR)/invtoj.o, $(wildcard $(CPDIR)/*.o))
CPLIBS = -L../local/lib -lzp_poly -lff_poly $(LIBS)

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
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

cm_method:  cm_method.o cornacchia.o fproot.o cminv.o invj.o class_inv_mpz.o
	$(CC) $(CFLAGS) -o $@ $^ $(CPOBJS) $(CPLIBS)

cmjtest:    cmjtest.o cminv.o invj.o fproot.o class_inv_mpz.o
	$(CC) $(CFLAGS) -o $@ $^ $(CPOBJS) $(CPLIBS)

zpbench:    zpbench.o fproot.o
	$(CC) $(CFLAGS) -o $@ $^ -L../local/lib -lzp_poly $(LIBS)

class_inv_mpz.o: $(CPDIR)/class_inv_mpz.c
	$(CC) $(CFLAGS) -c $< -o $@

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
