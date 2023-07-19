#!/bin/bash

g++ ~/code/phys/Crank-Nicolson2D-cxx/nicolson2D-rand.cxx -o nicolson2D-rand-cxx -I/usr/include/eigen3 -O3
~/code/phys/nicolson2D-rand-cxx
#python3 ~/code/phys/Crank-Nicolson2D-cxx/plot.py
