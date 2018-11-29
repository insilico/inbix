#!/bin/bash
#
# build_release.sh - Bill White - 11/29/18

cmake -DCMAKE_VERBOSE_MAKEFILE="OFF" -DCMAKE_BUILD_TYPE="Release" ..
make clean
make -j4

