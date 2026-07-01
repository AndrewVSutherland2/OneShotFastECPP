CC = gcc
CFLAGS = -O3 -march=native -m64 -funroll-loops -std=gnu11 -Wall -Wextra -fopenmp
LIBS = -lgmp -lm

PROGRAMS = cmsearch dscan smoothtest roottest hdroot

all: $(PROGRAMS)

cmsearch: cmsearch.o cornacchia.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

dscan: dscan.o cornacchia.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

smoothtest: smoothtest.o smooth.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

roottest: roottest.o fproot.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

hdroot: hdroot.o fproot.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

cornacchia.o cmsearch.o dscan.o: cornacchia.h
smooth.o smoothtest.o: smooth.h
fproot.o roottest.o hdroot.o: fproot.h

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o $(PROGRAMS)
