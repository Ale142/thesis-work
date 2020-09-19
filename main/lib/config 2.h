/* Define for simulation type goes here */ 
#define N_VERTEX 100
#define TIMER 1000
#define N_RAN 100
#define STEP

#ifdef STEP
/* Number of adjacents extracted if is performed STEP*/
#define K 1

#endif

// Probability that a ball will die.
#define DEATH_PROB 0.01

/* Type of the graph */
// #define FULL_GRAPH
// #define REGULAR_GRAPH
#define DEGREE_SEQUENCE

#if defined(REGULAR_GRAPH) || defined(DEGREE_SEQUENCE)
    #define DEGREE1 10
    #define DEGREE2 50
#endif