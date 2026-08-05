#!/bin/bash -e

source /opt/abacus-toolchain/install/setup
cmake -B build -G Ninja \
    -DBUILD_TESTING=ON \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
    -DENABLE_MLALGO=ON \
    -DENABLE_LIBXC=ON \
    -DENABLE_LIBRI=ON \
    -DENABLE_GOOGLEBENCH=ON \
    -DENABLE_RAPIDJSON=ON \
    -DENABLE_FLOAT_FFTW=ON \
    -DENABLE_DFTD4=ON \
    -DNEP_DIR=/opt/abacus-toolchain/install/NEP_CPU-main
cmake --build build -j $(nproc)
cmake --install build
