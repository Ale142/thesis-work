GSL = -I/usr/local/include 
IGRAPH = -I/usr/local/Cellar/igraph/0.7.1_6/include/igraph 
CC=gcc
SOURCES=simulazione.c
EXECUTABLE=simulazione
# CFLAGS= -c -Wall -Ic:/gsl
IGRAPH_LIB = -ligraph -L/usr/local/Cellar/igraph/0.7.1_6/lib
LIBS= -L/usr/local/lib -lgsl -lgslcblas -lm

run: bin/main
	bin/main

simulazione: build/simulazione.o build/frequency.o build/theoric.o
	build/simulazione.o

build/theoric.o: src/theoric.c lib/theoric.h
	gcc $(GSL) $(IGRAPH) -c $< -o $@

build/frequency.o: src/frequency.c lib/frequency.h
	gcc $(GSL) $(IGRAPH) -c $< -o $@

build/simulazione.o: src/simulazione.c lib/config.h
	gcc $(GSL) $(IGRAPH) -c $< -o $@

build/main.o: src/main.c lib/config.h
	gcc $(GSL) $(IGRAPH) -c $< -o $@

bin/main:  build/main.o build/simulazione.o build/frequency.o build/theoric.o
	gcc  -o bin/main build/theoric.o build/frequency.o build/simulazione.o build/main.o $(IGRAPH_LIB) $(LIBS)

clean:
	rm -f build/* bin/*

hard-clear:
	find . -type f -name '*.dat' -delete