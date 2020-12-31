#!/usr/bin/env bash

make clean
autoreconf -iv --install
./configure CC=gcc
make

make install
make dist
