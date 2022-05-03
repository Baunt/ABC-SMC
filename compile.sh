#!/bin/sh

#export CC=/usr/bin/clang-14
#export CXX=/usr/bin/clang++-14

rm -r build
mkdir -p build/
cd build
cmake ../
cmake --build .
cd ..
cp build/bin/ABC_SMC_Algorithm .
