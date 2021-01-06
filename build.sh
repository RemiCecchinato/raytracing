#!/bin/bash

DEBUG="-Og -g"

gcc src/main.c -o main -lm $DEBUG