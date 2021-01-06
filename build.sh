#!/bin/bash

DEBUG="-g3"

gcc src/main.c -o main -lm $DEBUG -O3 -lpthread -D_REENTRANT -DUSE_THREADS=1