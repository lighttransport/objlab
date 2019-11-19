#!/bin/bash


rm -rf build
mkdir build
cd build

CXX=clang++ CC=clang cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 ..

mv compile_commands.json ..
cd ..
