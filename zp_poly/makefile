# In-tree makefile for zp_poly (adapted from the original, which used absolute
# paths under /home/drew/dev).  Builds libzp_poly.a here and installs it plus the
# headers into the project-local prefix (default ../local), matching ff_poly.
#
#   make            # build the static library
#   make install    # stage into $(PREFIX)/lib and $(PREFIX)/include/zp_poly
#   make tests      # build the library's own test/benchmark programs
#   make clean

CC = gcc
CFLAGS = -O3 -march=native -m64 -fomit-frame-pointer -funroll-loops -std=gnu99 -Wall
INCLUDES = -I.
LIBS = -lgmp -lm
PREFIX ?= $(abspath $(CURDIR)/../local)

HEADERS = zp_poly_main.h zp_poly_mul.h zp_poly_util.h zp_poly_pair.h zp2.h zp_roots.h
SOURCES = zp_poly_small.c zp_poly_util.c zp_poly_mul.c zp_poly_div.c zp_poly_gcd.c \
          zp_poly_exp.c zp_poly_comp.c zp2.c zp_roots.c
OBJECTS = $(SOURCES:.c=.o)

# test/benchmark mains embedded in the sources behind -D flags
TESTS = zp_poly_find_split_root_test zp_poly_find_split_root_cost zp_poly_find_root_test \
        zp_poly_mul_test zp_poly_div_test zp_poly_gcd_test

all: libzp_poly.a

libzp_poly.a: $(OBJECTS)
	ar r $@ $(OBJECTS)
	ranlib $@

install: all
	mkdir -p $(PREFIX)/lib $(PREFIX)/include/zp_poly
	cp -v zp_poly.h $(PREFIX)/include
	cp -v $(HEADERS) $(PREFIX)/include/zp_poly
	cp -v libzp_poly.a $(PREFIX)/lib

tests: $(TESTS)

zp_poly_find_split_root_test: zp_poly_exp.c libzp_poly.a
	$(CC) $(CFLAGS) $(INCLUDES) -DZP_POLY_FIND_SPLIT_ROOT_TEST -o $@ $< libzp_poly.a $(LIBS)

zp_poly_find_split_root_cost: zp_poly_exp.c libzp_poly.a
	$(CC) $(CFLAGS) $(INCLUDES) -DZP_POLY_FIND_SPLIT_ROOT_COST -o $@ $< libzp_poly.a $(LIBS)

zp_poly_find_root_test: zp_poly_exp.c libzp_poly.a
	$(CC) $(CFLAGS) $(INCLUDES) -DZP_POLY_FIND_ROOT_TEST -o $@ $< libzp_poly.a $(LIBS)

zp_poly_mul_test: zp_poly_mul.c libzp_poly.a
	$(CC) $(CFLAGS) $(INCLUDES) -DZP_POLY_MUL_TEST -o $@ $< libzp_poly.a $(LIBS)

zp_poly_div_test: zp_poly_div.c libzp_poly.a
	$(CC) $(CFLAGS) $(INCLUDES) -DZP_POLY_DIV_TEST -o $@ $< libzp_poly.a $(LIBS)

zp_poly_gcd_test: zp_poly_gcd.c libzp_poly.a
	$(CC) $(CFLAGS) $(INCLUDES) -DZP_POLY_GCD_TEST -o $@ $< libzp_poly.a $(LIBS)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

clean:
	rm -f $(OBJECTS) libzp_poly.a $(TESTS)
