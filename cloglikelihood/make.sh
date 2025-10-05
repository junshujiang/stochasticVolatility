#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 <file.c>"
  exit 1
fi

SRC="$1"
OUTDIR="./"

if [ ! -f "$SRC" ]; then
  echo "No such file: $SRC"
  exit 1
fi

# INLA include 路径
INC="$(Rscript -e 'cat(system.file("include", package="INLA"))')"
if [ ! -f "$INC/cgeneric.h" ]; then
  echo "ERROR: cgeneric.h not found under: $INC"
  echo "Check INLA installation: install.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/stable'))"
  exit 1
fi

# 动态库后缀（Linux/mac = .so, Windows = .dll）
DYN_EXT="$(Rscript -e 'cat(.Platform$dynlib.ext)')"

mkdir -p "$OUTDIR"

BASENAME="$(basename "$SRC" .c)"
OBJ="$OUTDIR/$BASENAME.o"
SO="$OUTDIR/$BASENAME$DYN_EXT"

# 允许外部覆盖 CFLAGS/LDFLAGS
CFLAGS="${CFLAGS:- -Wall -O2 -fPIC}"
LDFLAGS="${LDFLAGS:- -shared -lm}"

# 编译
gcc -I"$INC" $CFLAGS -c "$SRC" -o "$OBJ"
# 链接
gcc $LDFLAGS -o "$SO" "$OBJ"
# 清理
rm -f "$OBJ"

echo "Built: $SO"