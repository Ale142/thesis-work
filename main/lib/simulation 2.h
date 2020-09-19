#ifndef __SIMULATION__
#define __SIMULATION__
#include <igraph/igraph.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <gsl/gsl_cdf.h>
#include <sys/time.h>


// #include "../lib/config.h"

/* Check if x is null */
#define NULL_CHECK(x) ({int val; if(x == 0) val = 1; else val = 0; val;})


typedef enum TYPE_RESULT{FREQUENCY, THEORIC} TYPE_RESULT;
typedef enum TYPE_OF_GRAPH{FULL_GRAPH, REGULAR_GRAPH, DEGREE_SEQUENCE, WATTS_STROGATZ, BARABASI} TYPE_OF_GRAPH;
typedef enum TYPE_OF_THEORIC {BINOMIAL, GEOMETRIC, GEOM_GEOM_1, WITH_CHOOSE, STABLE_DEGREE} TYPE_OF_THEORIC;

/* Global variables for simulation */
int *palline, *frequency,  *zero_balls;
double **frequency_nran_degree1, **frequency_nran_degree2;
char *file_name,  *directory, *new_file_name;
double **frequency_nran;

/* Global variable of graph */
igraph_t graph;

/* Struct for Watts Strogatz graph*/
extern struct config_watts_strogatz {
    int dimension;
    int neighborhood;
    double rewire_p;
    int loops;
    int multiple;
}config_watts_strogatz;
/* All parameters of simulation */
extern struct config_simulation {
    int timer;
    int t_trans;
    int n_vertex;
    int n_ran;
    int k;
    int with_step;
    int death_phase;
    float death_prob;
    char type_of_graph[128];
    TYPE_OF_GRAPH type;
    int degree1;
    int degree2;
    char *file_name_frequency;
    char *file_name_theoric;
    
}config_simulation;

struct config_simulation *conf;
struct config_watts_strogatz *conf_watts;

igraph_vector_t get_all_vertex(igraph_t graph);
igraph_vector_t get_adjacent(igraph_t graph, int vertex);
int number_of_vertices_graph(const igraph_t graph);
int get_degree(int vertex_id);
void print_status(int *palline);

void death_phase1(gsl_rng* random, igraph_vector_t all_vertex);
void death_phase2(gsl_rng* random, igraph_vector_t all_vertex);
void simulation(int seed,  igraph_t graph);
void multiple_simulation(igraph_t graph);
int find_min(int* neighbors_extracted, int vertex);

//for degree sequence
void inizialize_vector(igraph_vector_t vector, int from, int to, int value);

//array's stuff
void array_copy_int_to_double(int *array1, double *array2, int length);

void get_file_name(TYPE_RESULT T);
void write_parameters(FILE *output);

#endif