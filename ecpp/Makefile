CC = gcc
# class_inv.h (invariant codes + inv_string) and ff_poly.h are needed by cminv/cm_method
CPINC = -isystem ../classpoly_v1.0.3 -isystem ../local/include
CFLAGS = -O3 -march=native -m64 -funroll-loops -std=gnu11 -Wall -Wextra -fopenmp $(CPINC)
LIBS = -lgmp -lm

PROGRAMS = cmsearch dscan smoothtest roottest hdroot cm_method oneshot
TESTPROGS = invjtest cmjtest asmtest

all: $(PROGRAMS)
tests: $(TESTPROGS)

cmsearch:   cmsearch.o cornacchia.o
dscan:      dscan.o cornacchia.o
smoothtest: smoothtest.o smooth.o
roottest:   roottest.o fproot.o
hdroot:     hdroot.o fproot.o
cm_method:  cm_method.o cornacchia.o fproot.o cminv.o invj.o
oneshot:    oneshot.o cornacchia.o smooth.o fproot.o curve.o
invjtest:   invjtest.o invj.o fproot.o
cmjtest:    cmjtest.o cminv.o invj.o fproot.o
asmtest:    asmtest.o curve.o fproot.o cornacchia.o smooth.o

$(PROGRAMS) $(TESTPROGS):
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

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
