#include "../lib/simulation.h"
#include "../lib/theoric.h"

#define READ_INT "\t%*c\t%d"
#define READ_FLOAT "\t%*c\t%f"
#define READ_DOUBLE "\t%*c\t%lf"
#define READ_STRING "\t%*c\t%s"

void read_conf_file(struct config_simulation* conf) {
    FILE *file = fopen("simulazione.conf", "r");
    char variable[128];
    while (fscanf(file, "%s", variable) != EOF) {
        if(strcmp(variable, "TIMER") == 0) {
            fscanf(file, READ_INT, &conf->timer);
        } else if(strcmp(variable, "T_TRANS") == 0) {
            fscanf(file, READ_INT, &conf->t_trans);
        } else if(strcmp(variable, "N_VERTEX") == 0) {
            fscanf(file, READ_INT, &conf->n_vertex);
        } else if(strcmp(variable, "N_RAN") == 0) {
            fscanf(file, READ_INT, &conf->n_ran);
            frequency_nran = (double **) malloc(sizeof(double*) * conf->n_ran);

        } else if(strcmp(variable, "WITH_STEP") == 0) {
            fscanf(file, READ_INT, &conf->with_step);
        } else if(strcmp(variable, "K") == 0) {
            fscanf(file, READ_INT, &conf->k);
        } else if(strcmp(variable, "DEATH_PROB") == 0) {
            fscanf(file, READ_FLOAT, &conf->death_prob);
        } else if(strcmp(variable, "DEATH_PHASE") == 0) {
            fscanf(file, READ_INT, &conf->death_phase);
        } else if(strcmp(variable, "TYPE_OF_GRAPH") == 0) {
            fscanf(file, READ_STRING, conf->type_of_graph);
            if(strcmp(conf->type_of_graph, "FULL_GRAPH") == 0) 
                conf->type = FULL_GRAPH;
            else if(strcmp(conf->type_of_graph, "REGULAR_GRAPH") == 0)
                conf->type = REGULAR_GRAPH;
            else if(strcmp(conf->type_of_graph, "DEGREE_SEQUENCE") == 0) 
                conf->type = DEGREE_SEQUENCE;
            else if(strcmp(conf->type_of_graph, "WATTS_STROGATZ") == 0)
                conf->type = WATTS_STROGATZ;
            else if(strcmp(conf->type_of_graph, "BARABASI") == 0)
                conf->type = BARABASI;
        } else if(strcmp(variable, "DEGREE1") == 0) {
            fscanf(file, READ_INT, &conf->degree1);
        } else if(strcmp(variable, "DEGREE2") == 0) {
            fscanf(file, READ_INT, &conf->degree2);
        } 
    }

    fclose(file);
}


void print_config_simulation(struct config_simulation* conf) {
    printf("----------------------------------\n");
    printf("TIMER: %d\n", conf->timer);
    printf("TRANSITORY TIME: %d\n", conf->t_trans);
    printf("N_VERTEX %d\n", conf->n_vertex);
    if(NULL_CHECK(conf->n_ran) == 1) { printf("NULL N RAND");}
    
    printf("N_RAN %d\n", conf->n_ran);
    printf("K %d\n", conf->k);
    printf("DEATH_PROB %lf\n", conf->death_prob);
    printf("DEATH PHASE: %d\n", conf->death_phase);
    printf("TYPE OF GRAPH %s\n", conf->type_of_graph);
    if(strcmp(conf->type_of_graph, "DEGREE_SEQUENCE") == 0) {
        printf("DEGREE1 : %d\n", conf->degree1);
        printf("DEGREE2 : %d\n", conf->degree2);
    }
    printf("----------------------------------\n");

}

void read_config_watts_strogatz(struct config_watts_strogatz* conf) {
    FILE *config_file_graph = fopen("results/watts_strogatz/watts_strogatz.conf", "r");
    char variable[128];
     while (fscanf(config_file_graph, "%s", variable) != EOF) {
                if(strcmp(variable, "DIM") == 0) {
                    fscanf(config_file_graph, READ_INT, &conf->dimension);
                } else if (strcmp(variable, "NEI") == 0) {
                    fscanf(config_file_graph, READ_INT, &conf->neighborhood);
                } else if (strcmp(variable, "P") == 0) {
                    fscanf(config_file_graph, READ_DOUBLE, &conf->rewire_p);
                } else if (strcmp(variable, "LOOPS") == 0) {
                    fscanf(config_file_graph, READ_INT, &conf->loops);
                } else if (strcmp(variable, "MULTIPLE") == 0) {
                    fscanf(config_file_graph, READ_INT, &conf->multiple);
                }

            }
}

void generate_gnuplot_file() {
    FILE *gnuplot_file = fopen("run_simulation.p", "w");
    fprintf(gnuplot_file, 
        "reset\n"
        "set pointsize 1\n"
        "set grid xtics\n"
        "set grid ytics\n"
    );
    fprintf(gnuplot_file, "set title");
    write_parameters(gnuplot_file);

    // switch(conf->type) {
    //     case FULL_GRAPH:
    //     fprintf(gnuplot_file, "set title '# FULL GRAPH, conf->n_vertex : %d, N_RAN: %d, TIMER: %d T_TRANS: %d,,DEATH_PROB %f, DEATH_PHASE  %d, K %d' \n", conf->n_vertex, conf->n_ran, conf->timer,conf->t_trans, conf->death_prob, conf->death_phase, conf->k);
    //     break;
    //     case DEGREE_SEQUENCE:
    //     fprintf(gnuplot_file, "set title '# DEGREE SEQUENCE, conf->n_vertex : %d, N_RAN: %d, TIMER: %d,T_TRANS: %d, conf->degree1: %d, conf->degree2: %d, DEATH_PROB: %f, DEATH_PHASE: %d, K: %d '\n", conf->n_vertex, conf->n_ran, conf->timer,conf->t_trans, conf->degree1, conf->degree2, conf->death_prob, conf->death_phase, conf->k);
    //     break;
    //     case REGULAR_GRAPH:
    //         fprintf(gnuplot_file, "set title '# REGULAR GRAPH, conf->n_vertex : %d, N_RAN: %d, TIMER: %d,T_TRANS: %d, conf->degree1: %d, DEATH_PROB: %f, DEATH_PHASE: %d K: %d ' \n", conf->n_vertex, conf->n_ran, conf->timer,conf->t_trans, conf->degree1, conf->death_prob, conf->death_phase, conf->k);
    //     break;
    //     case WATTS_STROGATZ:
    //         fprintf(gnuplot_file, "set title '# WATTS-STROGATZ small-world, n_vertex: %d, N_RAN: %d, TIMER: %d,T_TRANS: %d, DEATH_PROB: %f, DEATH_PHASE: %d, K: %d;\\nDIMENSION LATTICE %d, NEIGHBOORHOOD: %d, REWIRE_P : %f, LOOPS: %d, MULTIPLE: %d '\n", conf->n_vertex, conf->n_ran, conf->timer,conf->t_trans, conf->death_prob, conf->death_phase, conf->k, conf_watts->dimension, conf_watts->neighborhood, conf_watts->rewire_p, conf_watts->loops, conf_watts->multiple);
    //     break;
    //     case BARABASI:
    //         fprintf(gnuplot_file, "set title '# BARABASI GRAPH, conf->n_vertex : %d, N_RAN: %d, TIMER: %d T_TRANS: %d,,DEATH_PROB %f, DEATH_PHASE  %d, K %d' \n", conf->n_vertex, conf->n_ran, conf->timer,conf->t_trans, conf->death_prob, conf->death_phase, conf->k);
    //     break;
    // }

    fprintf(gnuplot_file, "set xrange[:50]\n");
    fprintf(gnuplot_file, "plot '%s' u 1:3 with boxes t 'Frequency', '%s' u 1:2 with boxes t 'Theoric'\n", conf->file_name_frequency, conf->file_name_theoric);
    
}

int main(void) {

    igraph_vector_t out_deg;
    

    conf = (struct config_simulation*) malloc(sizeof(struct config_simulation) * 1);
    // zero_balls = (int *) malloc(sizeof(int) * (conf->timer + 1));
    read_conf_file(conf);

    printf("PARAMETERS of SIMULATION\n");
    print_config_simulation(conf);
    
    switch(conf->type) {
        case FULL_GRAPH:
            igraph_full(&graph, conf->n_vertex, IGRAPH_UNDIRECTED, 0);
        break;

        case REGULAR_GRAPH:
            /* Crea un grafo regolare: ogni vertice ha uguale grado */
            igraph_k_regular_game(&graph, conf->n_vertex, conf->degree1, IGRAPH_UNDIRECTED, IGRAPH_NO_MULTIPLE);
            // igraph_degree_sequence_game(&graph, )
        break;

        case DEGREE_SEQUENCE:


            igraph_vector_init(&out_deg, conf->n_vertex);
            inizialize_vector(out_deg, 0, conf->n_vertex / 2, conf->degree1);
            inizialize_vector(out_deg, conf->n_vertex / 2, conf->n_vertex, conf->degree2);


            igraph_degree_sequence_game(&graph, &out_deg, 0, IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE);
            
        break;

        case WATTS_STROGATZ: {
           

            char variable[128];
            conf_watts = (struct config_watts_strogatz*) malloc(sizeof(struct config_watts_strogatz) * 1);
           
            /* Read the config file*/
            read_config_watts_strogatz(conf_watts);

            printf("%d, %d, %lf, %d, %d\n", conf_watts->dimension, conf_watts->neighborhood, conf_watts->rewire_p, conf_watts->loops, conf_watts->multiple);

            igraph_watts_strogatz_game(&graph, conf_watts->dimension, conf->n_vertex, conf_watts->neighborhood, conf_watts->rewire_p, conf_watts->loops, conf_watts->multiple);
        }
        break;
        case BARABASI: {
            igraph_barabasi_game(&graph, 
            /* Numer of vertices */ conf->n_vertex,
            /* Power of preferential attachment*/ 1.0,
            /* m =  */ 2,
            /* outseq = */ 0,
            /* outpref = */ 0,
            /* A = */ 1.0,
            /* directed */ IGRAPH_UNDIRECTED,
            /* algorithm = */ IGRAPH_BARABASI_PSUMTREE,
            /* start_from = */ 0);

        }
    }

    printf("GRAPH GENERATED\n");
    printf("NUMBER OF VERTICES GRAPH %d\n", number_of_vertices_graph(graph));

    // printf("Start simulation...\n");
    printf("STARTING SIMULATION...\n");
    if(NULL_CHECK(conf->n_ran) == 0) { 
        /* Simulation is with multiple ran */
        multiple_simulation(graph);
    } else {
        /* Single Run simulation */
        simulation(0, graph);
    }

    /* Switch su fase di morte */
    switch(conf->death_phase) {
        case 1:
            // if(conf->type == DEGREE_SEQUENCE) {
            //     write_theoric_degree_seq_bin();
            // }
            ;
           // else
                switch (conf->with_step)
                {
                case 0:
                    /*WITH_STEP = 0*/
                    // write_theoric(BINOMIAL);
                    write_theoric(GEOM_GEOM_1);
                    break;
                
                case 1:
                    /* WITH_STEP = 1*/ 
                    write_theoric(WITH_CHOOSE);
                break;
               
                }
        break;
        case 2:
            /* Se la fase di morte è la seconda la distribuzione teorica è una geometrica:
                in caso il grafo fosse regolare/completo la geometrica è particolare: essa è
                una GEOM/GEOM/1
            */
            if(conf->with_step == 1) {
                if(conf->type == REGULAR_GRAPH || conf->type == FULL_GRAPH || conf->type == WATTS_STROGATZ || conf->type == BARABASI) {  
                    /* Se siamo nella fase in cui vengono scelti K vicini (WITH_STEP = 1) e nel caso il grafo sia con gradi 
                    omogenei, la teorica sarà calcolata in modo particola (v. documento pdf twochoices)*/  
                    write_theoric(WITH_CHOOSE);
                } else
                {
                    write_theoric(STABLE_DEGREE);
                    write_theoric(WITH_CHOOSE);
                }
                
            } else {
                /* Altrimenti controlliamo semplicemente il tipo di grafo*/
                if(conf->type == REGULAR_GRAPH || conf->type == FULL_GRAPH || conf->type == WATTS_STROGATZ || conf->type == BARABASI) {
                    /* Se il grafo ha gradi omogenei la teorica è una GEOM/GEOM/1*/
                    write_theoric(GEOM_GEOM_1);
                }
                else {
                    /* Qui il grafo è un DEGREE_SEQUENCE (ho dei nodi stabili e instabili)
                        Calcolo una geometrica con le pdf dei nodi stabili
                    */
                    write_theoric(STABLE_DEGREE);
                    write_theoric(GEOM_GEOM_1);
                    
                }
            }
        break;
    }

    
    
    generate_gnuplot_file();
    printf("GENERATED run_simulation.p; \nFor run simulation open gnuplot in the main directory and digit load 'run_simulation.p' \n");

    free(directory);
    free(file_name);

    printf("END SIMULATION\n");
}
