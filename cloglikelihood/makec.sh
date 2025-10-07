#!/usr/bin/env bash
set -euo pipefail

# 可选：传一个输出目录；不传则复制到当前目录


BUILD_DIR="build"
ROOT_DIR="$(pwd)"

# 干净重建
rm -rf "$BUILD_DIR"
mkdir -p "$BUILD_DIR"


echo "[1/4] CMake configure..."
cd "$BUILD_DIR"
cmake -DCMAKE_BUILD_TYPE=Release ..

echo "[2/4] Build..."
# 可按需把 -j 后的并行度改大/小
cmake --build . -j

echo "[3/4] Copy artifacts..."


mv -f "gaussiantest.so" ../

echo "[4/4] Cleanup build dir..."
cd "$ROOT_DIR"
rm -rf "$BUILD_DIR"

echo "Done."
