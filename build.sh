#!/bin/bash

DEBUG="-g3 -DUSE_THREADS=0"
RELEASE="-O3 -DUSE_THREADS=1"

OPT=$RELEASE

gcc src/main.c -o main -lm $OPT -lpthread -D_REENTRANT