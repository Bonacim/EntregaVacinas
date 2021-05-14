#!/bin/bash

# Points.
npts=$(head -1 input.txt); tail -$npts input.txt > input2.txt; gnuplot points.plot
rm input2.txt;

# Tree.
gnuplot tree.plot

# Cycle.
gnuplot cycle.plot

# Merge.
#pdfjam --nup 1x3 pontos.pdf arvore.pdf ciclo.pdf --outfile plots.pdf
