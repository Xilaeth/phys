#!/bin/bash

g++ ~/code/phys/Crank-Nicolson-k-space-cxx/nicolson.cxx -o nicolson-cxx -I/usr/include/eigen3 -O3
~/code/phys/nicolson-cxx
#python3 ~/code/phys/Crank-Nicolson-cxx/plot.py
