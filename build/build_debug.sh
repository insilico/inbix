#!/bin/bash
#
# build_debug.sh - Bill White - 8/5/18

cmake -DCMAKE_VERBOSE_MAKEFILE="ON" -DCMAKE_BUILD_TYPE="Debug" ..
make clean
make
