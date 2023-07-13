#!/bin/bash

g++ ~/code/phys/Crank-Nicolson2D-cxx/nicolson2D.cxx -o nicolson2D-cxx -I/usr/include/eigen3 -O3
~/code/phys/nicolson2D-cxx
python3 ~/code/phys/Crank-Nicolson2D-cxx/plot.py
