# Top-level build for the OneShotFastECPP dependencies.
#
# Everything builds and installs *within this project tree* -- nothing is
# written to /usr/local.  ff_poly is compiled and "installed" into ./local
# (a project-local prefix), and classpoly is then compiled against it.
#
#   make              # build ff_poly + classpoly + the ecpp tools (incl. oneshot)
#   make test         # build, then run ALL unit tests over |D| <= 1000
#   make test MAXD=200 # run the tests over a smaller range (faster)
#   make test-quick   # run the tests over |D| <= 200
#   make clean        # remove object files and binaries
#   make distclean    # also remove the ./local prefix and ./work scratch dirs

ROOT   := $(abspath $(CURDIR))
PREFIX := $(ROOT)/local
FFDIR  := $(ROOT)/ff_poly_v2.0.0
ZPDIR  := $(ROOT)/zp_poly
CPDIR  := $(ROOT)/classpoly_v1.0.3
WORK   := $(ROOT)/work
MAXD   ?= 1000

ECPPDIR := $(ROOT)/ecpp

.PHONY: all ff_poly zp_poly classpoly ecpp dirs test test-quick clean distclean

all: classpoly ecpp

# ---- ecpp: the ECPP tools (dscan, smoothtest, roottest, cm_method, oneshot) ----
# built against classpoly's headers/objects and the staged ff_poly + zp_poly.
ecpp: classpoly zp_poly
	$(MAKE) -C $(ECPPDIR)

# ---- ff_poly: build the static library and stage it under $(PREFIX) ----
ff_poly:
	$(MAKE) -C $(FFDIR)
	$(MAKE) -C $(FFDIR) install PREFIX=$(PREFIX)

# ---- zp_poly: large-p F_p[x] library (Harvey-Sutherland); used by classpoly's
# mpz_j_from_inv (class_inv_mpz.c), which the ecpp tools link ----
zp_poly:
	$(MAKE) -C $(ZPDIR)
	$(MAKE) -C $(ZPDIR) install PREFIX=$(PREFIX)

# ---- classpoly: build against the staged ff_poly ----
classpoly: ff_poly dirs
	$(MAKE) -C $(CPDIR) PREFIX=$(PREFIX)

# ---- default class-polynomial output directory (CRT scratch is per-process /tmp) ----
dirs:
	mkdir -p $(WORK)/H_files

# Run the full test suite (Test 1: classpoly vs PARI; Test 2: invariant->j checks).
test: classpoly
	cd $(ROOT)/tests && ./run_tests.sh $(MAXD)

test-quick: classpoly
	cd $(ROOT)/tests && ./run_tests.sh 200

clean:
	$(MAKE) -C $(FFDIR) clean
	$(MAKE) -C $(ZPDIR) clean
	$(MAKE) -C $(CPDIR) clean
	$(MAKE) -C $(ECPPDIR) clean

distclean: clean
	rm -rf $(PREFIX) $(WORK)
