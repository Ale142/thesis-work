/*
Mediare i risultati N RAN
Sommo i vari valori e divido per conf->n_ran

gsl_cdf_binomial_Q crea una funzione di ripartizione con valoi decrescenti
da 1
gsl_cdf_binomial_P crea una funzione di ripartizione con valori crescenti



TODO (26-03-20): ok
    Controllare se il bilanciamento è ancora più evidente:
         1. estraggo un nodo a caso (uniform int)
         2. estraggo un suo vicino a caso (uniform int)
         3. ripeto 2 fino ad un certo parametro K (eventualmente
            può capitare che estragga più volte gli stessi nodi)
         4. Trovo il nodo (tra quelli estratti nei passi 3-4)
            che ha il minor numero di palline e inserisco la pallina in quel nodo

    Palline che muoiono

    Ad ogni istante (ogni iterazione di conf->timer) una pallina muore con probabilità Q (parametro del sistema)
    1° fase, fase di morte delle palline:
        si esaminano tutti i nodi e guardo per ognuno il numero di palline che hanno; es: osservo che il nodo Ni
        ha M palline: genero una binomiale (gsl_ran_binomial) di parametri M e Q -> T = questo mi fornisce il numero
        di palline che muoiono per quel nodo Ni. Sottraggo M - T = nuovo numero di palline del nodo Ni.
    2° fase, fase di scelta del nodo e inserimento palline (uguale a prima)

    => Cosa succede: il sistema e il numero di palline non tente ad infinito ma a un numero stabile di palline


TODO (15-04-20)
    Devo suddividere in 2 array il conteggio delle frequency quando uso grado DEGREE_SEQUENCE: così poi
    posso facilmente calcolare la curva delle pdf solo di un certo grado(quello stabile).
    Per normalizzare le frequenze devo dividere per il N° di nodi che hanno quel grado analizzato ( e non per 
    N_VERTEX!!!). 
    La probabilità per calcolare la teorica è invece grado_vertice/somma di tutti i gradi


     
*/

#include "../lib/simulation.h"
#include "../lib/frequency.h"
// #define REGULAR_GRAPH
// #define conf->n_ran 10


/*
    Imposta il file_name e la destinazione in base ai parametri definiti nel file config.h
*/
void get_file_name(TYPE_RESULT T) {
    file_name = (char *)malloc(sizeof(char) * 100);
    directory = (char *)malloc(sizeof(char) * 100);

    switch(T) {
        case FREQUENCY:
            
            if(NULL_CHECK(conf->n_ran == 0)) {
                
                switch(conf->type) {
                case DEGREE_SEQUENCE:
                    sprintf(file_name, "avg_deg_seq_frequency_%d.dat", conf->timer);
                    sprintf(directory, "results/degree_seq/frequency/");
                break;
                case REGULAR_GRAPH:
                    sprintf(file_name, "avg_regular_graph_frequency_%d.dat", conf->timer);
                    sprintf(directory, "results/regular_graph/frequency/"); 
                break;
                case FULL_GRAPH:
                    sprintf(file_name, "avg_full_graph_frequency_%d.dat", conf->timer);
                    sprintf(directory, "results/full_graph/frequency/");
                break;
                case WATTS_STROGATZ:
                    sprintf(file_name, "avg_watts_strogats_frequency_%d.dat", conf->timer);
                    sprintf(directory, "results/watts_strogatz/frequency/");
                break;
                case BARABASI:
                    sprintf(file_name, "avg_barabasi_frequency_%d.dat", conf->timer);
                    sprintf(directory, "results/barabasi/frequency/");
                break;
                }
            } else {
                switch(conf->type) {
                    case DEGREE_SEQUENCE:
                        sprintf(file_name, "deg_seq_frequency_%d.dat", conf->timer);
                        sprintf(directory, "results/degree_seq/frequency/");
                    break;
                    case REGULAR_GRAPH: 
                        sprintf(file_name, "regular_graph_frequency_%d.dat", conf->timer);
                        sprintf(directory, "results/regular_graph/frequency/");
                    break;
                    case FULL_GRAPH:
                        sprintf(file_name, "full_graph_frequency_%d.dat", conf->timer);
                        sprintf(directory, "results/full_graph/frequency/");
                    break;
                    case WATTS_STROGATZ:
                        sprintf(file_name, "watts_strogats_frequency_%d.dat", conf->timer);
                        sprintf(directory, "results/watts_strogatz/frequency/");
                    break;
                    case BARABASI:
                        sprintf(file_name, "barabasi_frequency_%d.dat", conf->timer);
                        sprintf(directory, "results/barabasi/frequency/");
                    break;
                }
            }
            if(conf->with_step == 1) {
                new_file_name = (char *)malloc(sizeof(char) * 100);
                strcpy(new_file_name, "step_");
                strcat(new_file_name, file_name);
                strcpy(file_name, new_file_name);
            }
            if(conf->t_trans != 0) {
                new_file_name = (char *) malloc(sizeof(char) * 100);
                strcpy(new_file_name, "trans_");
                strcat(new_file_name, file_name);
                strcpy(file_name, new_file_name);
            }
        conf->file_name_frequency = (char*) malloc(sizeof(directory) * strlen(directory));
        break;
        case THEORIC:
                switch(conf->type) {
                    case DEGREE_SEQUENCE:
                        sprintf(file_name, "deg_seq_theoric_%d.dat", conf->timer);
                        sprintf(directory, "results/degree_seq/theoric/");
                    break;
                    case REGULAR_GRAPH:
                        sprintf(file_name, "regular_graph_theoric_%d.dat", conf->timer);
                        sprintf(directory, "results/regular_graph/theoric/");
                    break;
                    case FULL_GRAPH:
                        sprintf(file_name, "full_graph_theoric_%d.dat", conf->timer);
                        sprintf(directory, "results/full_graph/theoric/");
                    break;
                    case WATTS_STROGATZ:
                        sprintf(file_name, "watts_strogatz_%d.dat", conf->timer);
                        sprintf(directory, "results/watts_strogatz/theoric/");
                    break;
                    case BARABASI:
                        sprintf(file_name, "barabasi_%d.dat", conf->timer);
                        sprintf(directory, "results/barabasi/theoric/");
                    break;
            }
        conf->file_name_theoric = (char*) malloc(sizeof(directory) * strlen(directory));
        break;
    }

    /* Concat direcotry name and file name: the path's file to write results is stored in directory*/
    // strcat(directory, file_name);
}


/* Scrive su file i parametri della simulazione */
void write_parameters(FILE *output) {
    switch(conf->type) {
        case FULL_GRAPH:
        fprintf(output, "'# FULL GRAPH, conf->n_vertex : %d, N_RAN: %d, TIMER: %d, T_TRANS: %d ,DEATH_PROB %f, DEATH_PHASE  %d, K %d'\n", conf->n_vertex, conf->n_ran, conf->timer,conf->t_trans, conf->death_prob, conf->death_phase, conf->k);
        break;
        case DEGREE_SEQUENCE:
        fprintf(output, "'# DEGREE SEQUENCE, conf->n_vertex : %d, N_RAN: %d, TIMER: %d,T_TRANS: %d, conf->degree1: %d, conf->degree2: %d, DEATH_PROB: %f, DEATH_PHASE: %d, K: %d'\n", conf->n_vertex, conf->n_ran, conf->timer,conf->t_trans, conf->degree1, conf->degree2, conf->death_prob, conf->death_phase, conf->k);
        break;
        case REGULAR_GRAPH:
            fprintf(output, "'# REGULAR GRAPH, conf->n_vertex : %d, N_RAN: %d, TIMER: %d,T_TRANS: %d, conf->degree1: %d, DEATH_PROB: %f, DEATH_PHASE: %d K: %d'\n", conf->n_vertex, conf->n_ran, conf->timer,conf->t_trans, conf->degree1, conf->death_prob, conf->death_phase, conf->k);
        break;
        case WATTS_STROGATZ:
            fprintf(output, "'# WATTS-STROGATZ small-world, n_vertex: %d, N_RAN: %d, TIMER: %d,T_TRANS: %d, DEATH_PROB: %f, DEATH_PHASE: %d, K: %d; DIMENSION LATTICE %d, NEIGHBOORHOOD: %d, REWIRE_P : %f, LOOPS: %d, MULTIPLE: %d'\n", conf->n_vertex, conf->n_ran, conf->timer,conf->t_trans, conf->death_prob, conf->death_phase, conf->k, conf_watts->dimension, conf_watts->neighborhood, conf_watts->rewire_p, conf_watts->loops, conf_watts->multiple);
        break;
        case BARABASI:
            fprintf(output, "'# BARABASI, conf->n_vertex : %d, N_RAN: %d, TIMER: %d, T_TRANS: %d ,DEATH_PROB %f, DEATH_PHASE  %d, K %d'\n", conf->n_vertex, conf->n_ran, conf->timer,conf->t_trans, conf->death_prob, conf->death_phase, conf->k);
        break;
    }
}

/* Ritorna un vettore contenente tutti i vertici */
igraph_vector_t get_all_vertex(igraph_t graph) {
    igraph_vs_t all_vertex_selector;
    igraph_vit_t vertex_iterator;
    igraph_vector_t vertex_vector;

    int i = 0;
    // Seleziono tutti i vertici
    igraph_vs_all(&all_vertex_selector);

    // Creo iteratore su tutti i vertici
    igraph_vit_create(&graph, all_vertex_selector, &vertex_iterator);
    igraph_vector_init(&vertex_vector, IGRAPH_VIT_SIZE(vertex_iterator));
    while(!IGRAPH_VIT_END(vertex_iterator)) {
        igraph_vector_set(&vertex_vector, i, IGRAPH_VIT_GET(vertex_iterator));
        // vertex_vector[i] = IGRAPH_VIT_GET(vertex_iterator);
        IGRAPH_VIT_NEXT(vertex_iterator);
        i++;
    }
    return vertex_vector;
}

int number_of_vertices_graph(const igraph_t graph) {
    igraph_vs_t vs;
    igraph_integer_t *result;

    igraph_vs_all(&vs);
    igraph_vs_size(&graph,&vs, result);
    return *result;
}

/* Ritorna un vettore di tutti i vertici adiacenti al vertex passato come input*/
igraph_vector_t get_adjacent(igraph_t graph, int vertex) {
    igraph_vs_t adjacent_selector;
    igraph_vit_t adjacent_iterator;
    igraph_vector_t neighbors;
    int i = 0;

    //Seleziono adiacenti
    igraph_vs_adj(&adjacent_selector, vertex, IGRAPH_ALL);

    // Creo iteratore su adiacenti
    igraph_vit_create(&graph, adjacent_selector, &adjacent_iterator);
    igraph_vector_init(&neighbors, IGRAPH_VIT_SIZE(adjacent_iterator));
    while(!IGRAPH_VIT_END(adjacent_iterator)) {
        igraph_vector_set(&neighbors, i, IGRAPH_VIT_GET(adjacent_iterator));
        i++;
        IGRAPH_VIT_NEXT(adjacent_iterator);
    }
    return neighbors;
}

/* Ritorna il grado di un vertice */
int get_degree(int vertex_id) {
    igraph_vs_t adjacent_selector;
    igraph_integer_t degree;
    igraph_vs_adj(&adjacent_selector, vertex_id, IGRAPH_ALL);

    
    igraph_vs_size(&graph, &adjacent_selector, &degree);
    return degree;
}

void print_status(int *palline) {
    printf("STATUS URNE (vertex, number balls): \n");
    for(int i = 0;i < conf->n_vertex; i++) {
        printf("(%d , %d) \t", i, palline[i]);
    }
    printf("\n");
}

/* Fase 1 di morte: 
    Per ogni nodo del grafo estrggo un numero con probabilità binomiale
    (Bin (N_vertex, Death_prb)) che rappresenta quante palline devono morire per quel 
    nodo. Sottraggo quel numero al numero corrente di palline

 */
void death_phase1(gsl_rng* random, igraph_vector_t all_vertex) {
    int current_vertex, current_ball_of_vertex, number_of_balls_died, index;
    for(int i = 0; i < igraph_vector_size(&all_vertex); i++) {
        current_vertex = VECTOR(all_vertex)[i];
        current_ball_of_vertex = palline[current_vertex];

        // Calcolo il numero di palline che deve morire
        number_of_balls_died = gsl_ran_binomial(random, conf->death_prob, current_ball_of_vertex);

        // Aggiorno il nuovo numero di palline togliendo quelle morte
        palline[current_vertex] -= number_of_balls_died;
    }

    // printf("index %d\n", index);
    // current_vertex = VECTOR(all_vertex)[index_vertex_selected];
    // current_ball_of_vertex = palline[current_vertex];
    // number_of_balls_died = gsl_ran_binomial(random, conf->death_prob, current_ball_of_vertex);
    // palline[current_ball_of_vertex] -= number_of_balls_died;
}


/* Per ogni nodo del grafo, estraggo un numero casuale tra [0,1):
    se esso è minore della probabilità di morte allora tolgo una pallina
    dal nodo
*/
void death_phase2(gsl_rng* random, igraph_vector_t all_vertex) {
    int current_vertex;
    double p;
    for(int i = 0; i < igraph_vector_size(&all_vertex); i++) {
        current_vertex = VECTOR(all_vertex)[i];

        // Resituisco un numero compreso tra [0,1)
        p = gsl_rng_uniform(random);
        if(p < conf->death_prob && palline[current_vertex] > 0) 
            palline[current_vertex] -= 1;
        // Aggiorno il nuovo numero di palline togliendo quelle morte
       
    }
}

void array_copy_int_to_double(int *array1, double *array2, int len) {
    for(int i = 0; i < len; i++) {
        array2[i] = array1[i];
    }
}

void simulation(int seed,  igraph_t graph) {
    igraph_vector_t neighbors;
    igraph_vector_t all_vertex;
    struct timeval tv;

    gettimeofday(&tv,0);

    gsl_rng *random;


    // int time = conf->timer, *palline, *frequency;
    int time = conf->timer, timer_index = 0;
    long int my_seed = tv.tv_sec + tv.tv_usec;

    frequency = (int *) calloc(conf->timer + 1, sizeof(int));
    palline = (int *) calloc(conf->n_vertex, sizeof(int));
    zero_balls = (int *) calloc(conf->timer + 1, sizeof(int));

    random = gsl_rng_alloc(gsl_rng_rand);
    gsl_rng_set(random, my_seed);

    if(NULL_CHECK(conf->n_ran) != 0) {
        gsl_rng_set(random, seed + my_seed);
    }
    
    /* Seleziono tutti i vertici */
    all_vertex = get_all_vertex(graph);
    // igraph_vector_print(&all_vertex);

    // simulation
    while (timer_index < time) {

        printf("%d\n", timer_index);
        if(NULL_CHECK(conf->death_prob) != 1) {
        /* Simulation is with Death Probability*/
        // int index_vertex_selected_to_kill_balls = gsl_rng_uniform_int(random, igraph_vector_size(&all_vertex));
        // int selected_vertex_to_kill_balls = VECTOR(all_vertex)[index_vertex_selected_to_kill_balls];
        // int balls_to_kill = gsl_ran_binomial(random, conf->death_prob, palline[selected_vertex_to_kill_balls]);
        // palline[selected_vertex_to_kill_balls] -= balls_to_kill;
        // death_phase(random, all_vertex);
            switch(conf->death_phase) {
                case 1: 
                /* K palline vengono tolte da M tramite una distribuzione binomiale*/
                    death_phase1(random, all_vertex);
                break;
                case 2:
                /* Calcolo un P compreso tra [0,1): se questo è minore di DEATH_PROB tolgo UNA pallina al vertice*/
                    death_phase2(random, all_vertex);
                break;
            }
        }

        /* Estraiamo a caso uno dei vertici */
        int index_vertex_extracted = gsl_rng_uniform_int(random, igraph_vector_size(&all_vertex));
        int vertex_extracted = VECTOR(all_vertex)[index_vertex_extracted];
        // printf("INDEX EXTRACTED %d\t VERTEX ID %d\n", index_vertex_extracted, vertex_extracted);

        /*Selezione degli adiacenti del vertice selezionato prima*/
        neighbors = get_adjacent(graph, vertex_extracted);
        // printf("NEIGHBORS OF %d :\n", vertex_extracted);
        // igraph_vector_print(&neighbors);

        if(conf->with_step == 1) {
            /* Choose the vertex with mininum number of ball from conf->k neighbors */
            int *adjacents_extracted = (int *)malloc(sizeof(int) * conf->k);
            for(int i = 0; i < conf->k; i++) {

                /* Estraiamo un adiacente a caso e inseriamo la pallina in esso */
                int index_adjacent_extracted = gsl_rng_uniform_int(random, igraph_vector_size(&neighbors));
                // int adjacent_extracted = find_min(neighbors, vertex_extracted);
                adjacents_extracted[i] = VECTOR(neighbors)[index_adjacent_extracted];

                // if(palline[vertex_extracted] > palline[adjacent_extracted])
                //     palline[adjacent_extracted]++;
                // else
                //     palline[vertex_extracted]++;
            }
            int vertex_to_insert_ball = find_min(adjacents_extracted, vertex_extracted);
            palline[vertex_to_insert_ball]++;

        } else {
            // printf("here\n");
            // printf("EXTRACTED %d\n", adjacent_extracted);
            int index_adjacent_extracted = gsl_rng_uniform_int(random, igraph_vector_size(&neighbors));
            int adjcent_extracted = VECTOR(neighbors)[index_adjacent_extracted];
            // printf("INDEX ADJIACENT EXTRACTED %d\t ADJACENT ID %d\n", index_adjacent_extracted, adjcent_extracted);
            palline[adjcent_extracted]++;
            // printf("Inserted ball in vertex %d\n", adjcent_extracted);

        }
        count_k_balls(timer_index);
        
        if(timer_index > conf->t_trans) {
            count_frequency(frequency, palline);
            // printf("%d:\n ", time);
            //  print_frequency(frequency);
            //  printf("\n");
        }
        
       
        // if(conf->t_trans != 0)
        //     count_frequency_transitory();
        // else {
            
        // }
        // printf("---------------\n");
        // printf("\n");
        timer_index++;
    }

    if(NULL_CHECK(conf->n_ran) != 0) {
         // simulation(graph, neighbors, &palline, time, all_vertex);
        count_frequency(frequency, palline);

        // print_status(palline);
        // print_frequency(frequency);
        printf("WRITING RESULTS...\n");
        write_frequency(frequency);    
    }
}

/* Ritorna il vertice nell'array neighbors_extracted con il minor numero di palline */
int find_min(int* neighbors_extracted, int vertex) {
    int vertex_extracted = vertex, count = 0;
    // int *min = (int *) malloc(sizeof(int) * 0);
    for(int i = 0; i < conf->k; i++) {
        int neighbors = neighbors_extracted[i];
        if(palline[neighbors] < palline[vertex]) {
            vertex_extracted = neighbors;
        }
    }

    return vertex_extracted;
}


void multiple_simulation(igraph_t graph) {

    if(conf->type == DEGREE_SEQUENCE) {
        // Se il tipo di grafo è DEGREE_SEQUENCE, divido le frequenze in 2 array
        frequency_nran_degree1 = (double **) malloc(sizeof(double*) * conf->n_ran);
        frequency_nran_degree2 = (double **) malloc(sizeof(double*) * conf->n_ran);
    }

    int i = 0;
    while (i < conf->n_ran) {
        frequency_nran[i] = (double *) malloc(sizeof(double) * (conf->timer + 1));
        simulation(i + 1, graph);
        // count_frequency(frequency_nran[i], palline);
        
        // array_copy_int_to_double(frequency, frequency_nran[i], conf->timer + 1);
        // printf("after copy;");
        // for(int j = 0; j <= conf->timer; j++)
        //     printf("%f\t",frequency_nran[i][j]);
        // printf("RUN %d\n", i);
        // print_frequency(frequency);
        for(int j = 0; j <= conf->timer; j++) {
            frequency_nran[i][j] = (double) frequency[j] / (conf->timer - conf->t_trans);
        }

        // for(int j = 0; j <= conf->timer; j++) {
        //     printf("(%d - %f) \t", j, frequency_nran[i][j]);
        // }
        // printf("\n");
        // printf("\nAfter divide:");
        // for(int j = 0; j <= conf->timer; j++)
        //     printf("%f\t",frequency_nran[i][j]);
        // printf("\n");
       
        
        if(conf->type == DEGREE_SEQUENCE) {
            /* Se il tipo di grafo è un DEGREE_SEQUENCE calcolo le frequenze per i gradi 1 e 2*/
            frequency_nran_degree1[i] = (double *) malloc(sizeof(double) * (conf->timer + 1));
            frequency_nran_degree2[i] = (double *) malloc(sizeof(double) * (conf->timer + 1));
            
            count_frequency_single_degree(frequency_nran_degree1[i], palline, conf->degree1);
            count_frequency_single_degree(frequency_nran_degree2[i], palline, conf->degree2);
        }
        i++;
    }
        printf("END %d SIMULATION\n", conf->n_ran);
        printf("CALCULATE RESULTS...\n");
        double *average_frequency;
        average_frequency = (double *) malloc(sizeof(double) * (conf->timer + 1));

        // printf("RESULT FOR EACH RUN:\n");
        // for(int i = 0; i < conf->n_ran; i++){
        //     printf("RUN %d\n", i);
        //     for(int j = 0; j <= conf->timer; j++) {
        //         printf("%f\n", frequency_nran[i][j]);
        //     }
        //     printf("\n");
        // }

        calculate_average_frequency(average_frequency, frequency_nran);

        if(conf->type == DEGREE_SEQUENCE) {
            double *average_frequency_degree1, *average_frequency_degree2;
            average_frequency_degree1 = (double *) malloc(sizeof(double) * (conf->timer + 1));
            average_frequency_degree2 = (double *) malloc(sizeof(double) * (conf->timer + 1));
            calculate_average_frequency(average_frequency_degree1, frequency_nran_degree1);
            calculate_average_frequency(average_frequency_degree2, frequency_nran_degree2);

            write_single_frequency(average_frequency_degree1, average_frequency_degree2);

        }


        
        // printf("AVERAGE \n");
        // for(int h = 0; h < conf->timer; h++) {
        //     printf("(%d,%f)\t", h, average_frequency[h]);
        // }
        printf("WRITING AVERAGE FREQUENCY RESULTS...\n");
        write_average_frequency(average_frequency);

        printf("WRITE ZERO_BALLS DISTRIBUTIONS\n");
        write_k_balls();
        
}


void inizialize_vector(igraph_vector_t vector, int from, int to, int value) {
    for(int i = from; i < to; i++) {
        VECTOR(vector)[i] = value;
    }
}



