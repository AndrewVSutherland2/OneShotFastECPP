# Source this file to run classpoly with its data and output inside the project:
#
#     . ./setenv.sh
#     classpoly -23 0
#
# It points classpoly at the project's modular-polynomial files and a project
# output directory, and puts the classpoly and ecpp binaries (oneshot, cm_method,
# dscan, ...) on your PATH so they can find each other.
#
# CRT scratch files are handled by classpoly itself: it creates a private
# per-process directory under /tmp (override the base with CLASSPOLY_TMPDIR),
# so nothing needs to be set up here and concurrent runs never collide.

_oneshot_here="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"

export CLASSPOLY_PHI_DIR="$_oneshot_here/phi_files"   # modular polynomial files phi_*_*.txt
export CLASSPOLY_H_DIR="$_oneshot_here/work/H_files"  # default class-polynomial output dir
export ONESHOT_PCACHE_DIR="$_oneshot_here/work/pcache" # prime-product caches (survive reboots)
export PATH="$_oneshot_here/classpoly_v1.0.3:$_oneshot_here/ecpp:$PATH"

mkdir -p "$CLASSPOLY_H_DIR" "$ONESHOT_PCACHE_DIR"

echo "classpoly environment configured:"
echo "  CLASSPOLY_PHI_DIR = $CLASSPOLY_PHI_DIR"
echo "  CLASSPOLY_H_DIR   = $CLASSPOLY_H_DIR"
echo "  CRT scratch       = \$CLASSPOLY_TMPDIR or /tmp (per-process, auto-cleaned)"

unset _oneshot_here
