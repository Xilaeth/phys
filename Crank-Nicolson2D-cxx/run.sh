#!/bin/bash

g++ ~/code/phys/phys/Crank-Nicolson2D-cxx/main.cxx -o main -I/usr/include/eigen3
~/code/phys/phys/main
python ~/code/phys/phys/Crank-Nicolson2D-cxx/plot.py
