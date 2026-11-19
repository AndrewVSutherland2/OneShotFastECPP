#!/bin/bash
# Build every engine under test.  Run from the bundle directory on the box.
R=$(cd "$(dirname "$0")" && pwd)
cd "$R" || exit 1
mkdir -p src inst out
J=48
{
echo "== host"; hostname; cat /etc/instance-identity 2>/dev/null; lscpu | grep -E "Model name|^CPU\(s\)|Thread|Flags" | cut -c1-400; gcc --version | head -1; python3 --version

build_gmp () {   # name, extra configure args
    local name=$1; shift
    mkdir -p src/$name && tar xf gmp-6.3.0.tar.xz -C src/$name
    ( cd src/$name/gmp-6.3.0 && ./configure --prefix="$R/inst/$name" --disable-shared -q "$@" > "$R/out/$name.conf.log" 2>&1 \
      && make -j$J > "$R/out/$name.make.log" 2>&1 && make install > /dev/null 2>&1
      echo "$name: $(grep -m1 'path=' config.log)" )
}
build_gmp gmp-k8 &
build_gmp gmp-zen --build=zen3-pc-linux-gnu &
wait
ls -la inst/gmp-k8/lib/libgmp.a inst/gmp-zen/lib/libgmp.a

build_ecm () {   # name, gmp name
    local name=$1 gmp=$2
    mkdir -p src/$name && tar xzf ecm-7.0.6.tar.gz -C src/$name
    ( cd src/$name/ecm-7.0.6 && ./configure --with-gmp="$R/inst/$gmp" --disable-shared -q > "$R/out/$name.conf.log" 2>&1 \
      && make -j$J > "$R/out/$name.make.log" 2>&1
      ls -la ecm && cp ecm "$R/inst/$name" )
}
build_ecm ecm-k8 gmp-k8 &
build_ecm ecm-zen gmp-zen &
build_ecm ecm-zentuned gmp-zen &
wait
# the campaign binary: Ubuntu's gmp-ecm 7.0.5 against the system libgmp
mkdir -p inst/deb && for d in *.deb; do dpkg -x "$d" inst/deb; done
LD_LIBRARY_PATH=$R/inst/deb/usr/lib/x86_64-linux-gnu inst/deb/usr/bin/ecm -printconfig | head -4

build_avx () {   # name, march
    local name=$1 march=$2
    rsync -a avx-ecm-src/ src/$name/
    ( cd src/$name && make COMPILER=gcc ICELAKE=1 \
        OPT_FLAGS="-O3 -mavx -march=$march -DIFMA -DICELAKE -fcommon" \
        INC="-I. -I$R/inst/gmp-zen/include" LIBS="$R/inst/gmp-zen/lib/libgmp.a -lm -lpthread" \
        -j$J > "$R/out/$name.make.log" 2>&1
      ls -la avx-ecm && cp avx-ecm "$R/inst/$name" )
}
build_avx avx-icl icelake-client &
build_avx avx-native native &
wait
# tuned variant: let gmp-ecm measure its own thresholds on this host
( cd src/ecm-zentuned/ecm-7.0.6 && timeout 1500 make ecm-params > "$R/out/ecm-zentuned.tune.log" 2>&1
  echo "tune exit $?"; cat ecm-params.h 2>/dev/null | grep -v '^$' | head -40
  make -j$J > "$R/out/ecm-zentuned.make2.log" 2>&1; cp ecm "$R/inst/ecm-zentuned" )


for g in k8 zen; do gcc -O2 -I inst/gmp-$g/include -o inst/mb-$g mb.c inst/gmp-$g/lib/libgmp.a; done
gcc -O2 -o inst/mb-sys mb.c -lgmp 2>&1 | head -3
ls -la inst
echo BUILD_DONE
} > out/build.log 2>&1
