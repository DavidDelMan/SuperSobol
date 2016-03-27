#!/bin/bash

# g++ -O2 -std=c++0x SobolIndices.cpp SobolIndicesDriver.cpp Halton.cpp MT64.cpp InverseTransformation.cpp 

g++ -O2 -std=c++0x SuperSobolIndices.cpp SobolIndices.cpp SuperSobolDriver.cpp Halton.cpp MT64.cpp InverseTransformation.cpp MersenneTwister.cpp pdflib.cpp rnglib.cpp

# ./a.out 20000
# ./a.out 50000
# ./a.out 100000
# ./a.out 200000
# ./a.out 300000
# ./a.out 400000
# ./a.out 500000
# ./a.out 600000
# ./a.out 700000
# ./a.out 800000
# ./a.out 900000
# ./a.out 1000000

# ./a.out 25
# ./a.out 27.5
# ./a.out 30
# ./a.out 32.5
# ./a.out 35