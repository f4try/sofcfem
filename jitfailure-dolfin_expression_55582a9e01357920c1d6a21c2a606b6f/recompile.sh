#!/bin/bash
# Execute this file to recompile locally
c++ -Wall -shared -fPIC -std=c++11 -O3 -fno-math-errno -fno-trapping-math -ffinite-math-only -I/Users/zon/opt/anaconda3/envs/fenicsproject/include -I/Users/zon/opt/anaconda3/envs/fenicsproject/include/eigen3 -I/Users/zon/opt/anaconda3/envs/fenicsproject/.cache/dijitso/include dolfin_expression_55582a9e01357920c1d6a21c2a606b6f.cpp -L/Users/zon/opt/anaconda3/envs/fenicsproject/lib -L/Users/zon/opt/anaconda3/envs/fenicsproject/Users/zon/opt/anaconda3/envs/fenicsproject/lib -L/Users/zon/opt/anaconda3/envs/fenicsproject/.cache/dijitso/lib -Wl,-rpath,/Users/zon/opt/anaconda3/envs/fenicsproject/.cache/dijitso/lib -lpmpi -lmpi -lmpicxx -lpetsc -lslepc -lm -ldl -lz -lpthread -lhdf5 -lboost_timer -ldolfin -Wl,-install_name,/Users/zon/opt/anaconda3/envs/fenicsproject/.cache/dijitso/lib/libdijitso-dolfin_expression_55582a9e01357920c1d6a21c2a606b6f.so -olibdijitso-dolfin_expression_55582a9e01357920c1d6a21c2a606b6f.so