#ifndef __THEORIC__
#define __THEORIC__
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#pragma once 
#include "../lib/simulation.h"

void write_theoric(TYPE_OF_THEORIC type);
void write_theoric_bin(double *, double *);
void write_theoric_geom(double *, double *);
void write_theoric_degree_seq_bin();
void write_theoric_geom_geom_1(double *, double *);
#endif