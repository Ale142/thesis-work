#!/bin/zsh
gsl="-I/usr/local/include -L/usr/local/lib"
igraph="-I/usr/local/Cellar/igraph/0.7.1_6/include/igraph -L/usr/local/Cellar/igraph/0.7.1_6/lib"
# gcc $1 ${igraph} -ligraph -o run
gcc -Wall  $1 ${gsl} -lgsl -lgslcblas -lm ${igraph} -ligraph -o run