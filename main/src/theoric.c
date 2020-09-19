#include "../lib/theoric.h"

/* Curva teorica binomiale */
void write_theoric_bin(double *vector_pdf, double *vector_cdf) {
    for(int i = 0;i <= conf->timer; i++) {
        vector_pdf[i] = gsl_ran_binomial_pdf(i, 1.0 / conf->n_vertex ,conf->timer);
        vector_cdf[i] = gsl_cdf_binomial_P(i, 1.0 / conf->n_vertex, conf->timer);
    }
}

/* Curva teorica geometrica*/
void write_theoric_geo(double *vector_pdf, double *vector_cdf) {
    for(int i = 0;i <= conf->timer; i++) {
        vector_pdf[i] = gsl_ran_geometric_pdf(i, conf->death_prob);
        vector_cdf[i] = gsl_cdf_geometric_P(i, conf->death_prob);
    }
} 

/* Particolare tipo di geometrica (in caso di DEATH_PHASE 2)*/
void write_theoric_geom_geom_1(double *vector_pdf, double *vector_cdf) {
   
    /*
    La probabilità che un nodo abbia i nodi : 
        Pi(i) = (U/D)^ i - 1 * Pi(1),
        dove: U = probabilità di vita = p (1- q)
              D = probabilità di morte = q (1- p)
              Pi(1) = [p(q - p)] / [q^2 (1-p)]
              Pi(0) = [(q(1 - p)) / p] * Pi(1)
    p = 1/N_VERTEX, q = DEATH_PROB

    Per ragioni di stabilità la probabilità di vita(U) deve essere < probabilità di morte (D)
    */
    
   double U,D,p,q, Pi_0, Pi_1;
   p = 1.0 / conf->n_vertex;
   q = conf->death_prob;
   U = p * (1 - q);
   D = q * (1 - p);
   printf("P : %f, q: %f, U: %f, D:%f\n", p,q, U, D);
   printf("U : %f, D: %f\n", U, D);
   Pi_1 = ((p * (q - p)) / (pow(q,2.0) * (1 - p)));
   Pi_0 = (q * (1 - p) / p) * Pi_1;

   vector_pdf[0] = Pi_0 * Pi_1;
   for(int i = 0; i <= conf->timer; i++) {
       vector_pdf[i] = pow((U / D), (i-1)) * Pi_1;
       vector_cdf[i] = 0.0;
   }
   
}

/* Quando DEATH_PHASE = 2 e DEGREE_SEQUENCE
   Curva teorica per nodi stabili: serve per vedere che 
   alcuni nodi, con grado particolare, nel caso di DEGREE_SEQUENCE sono stabili

 */
void  write_theoric_stable_nodes(double *vector_pdf, double *vector_cdf) {
    /* Scrive la curva teorica e frequenza dei nodi stabili (quelli con grado minore) nel caso di grafo
       DEGREE_SEQUENCE
    */

   double N1 =  0.5 * conf->n_vertex; // Nodi di grado DEGREE 1
   double N2 =  0.5 * conf->n_vertex; // Nodi di grado DEGREE 2

   // Conto le frequenze sono per i primi conf->n_vertex / 2 nodi (quelli con grado DEGREE1, assumo che siano quelli stabili);
   double Q1 = (conf->degree1) / (N1 * conf->degree1 + N2 * conf->degree2); 
   double U,D,p,q, Pi_0, Pi_1;
   
   p = Q1;
   q = conf->death_prob;
   U = p * (1 - q);
   D = q * (1 - p);
   printf("P : %f, q: %f, U: %f, D:%f\n", p,q, U, D);
   printf("U : %f, D: %f\n", U, D);
   Pi_1 = ((p * (q - p)) / (pow(q,2.0) * (1 - p)));
   Pi_0 = (q * (1 - p) / p) * Pi_1;

   vector_pdf[0] = Pi_0 * Pi_1;
   for(int i = 0; i <= conf->timer; i++) {
       vector_pdf[i] = pow((U / D), (i-1)) * Pi_1;
       vector_cdf[i] = 0.0;
   }
}


/*
        Capitiamo in questo caso quando:
            DEATH_PHASE = 2
            WITH_ STEP = 1
            K = n
        ed il tipo di grafo è REGULAR_GRAPH oppure FULL_GRAPH

        (with_choose.pdf paper)
        
*/
void write_theoric_with_choose(double *vector_pdf, double *vector_cdf) {
    
   double U,D,p,q;
   int d = conf->k+1;
   
   p = 1.0 / conf->n_vertex; // Probilità di inserire palline
   
   q = conf->death_prob; // probabilità di morte
   U = p * (1 - q); // Calcolo p di vita (secondo slide)
   D = q * (1 - p); // Calcolo p di morte (secondo slide)
   printf("P : %f, q: %f, U: %f, D:%f\n", p,q, U, D);
   printf("U : %f, D: %f\n", U, D);
   for(int i = 0; i <= conf->timer; i++) {
       /* s_i = (U/D)^ ((d^i-1)/(d-1)) */
       double d_i = pow(d, i) - 1;
       d_i = d_i / (d-1);
       vector_cdf[i] = pow(U/D, (d_i));
   }
   
   for(int i = 0; i < conf->timer; i++) {
       vector_pdf[i] = vector_cdf[i] - vector_cdf[i+1];
   }
   vector_pdf[conf->timer] = vector_cdf[conf->timer];
}

/* Curva teorica per degree_sequence graph*/
void write_theoric_degree_seq_bin() {
    double *vector_pdf = (double *) malloc(sizeof(double) * (conf->timer + 1));
    // double *vector_cdf = (double *) malloc(sizeof(double) * conf->timer);
    get_file_name(THEORIC);
    strcat(directory, file_name);
    FILE *output = fopen(directory, "w");
    write_parameters(output);
    fprintf(output, "# %s %d, %s %d\n", "DEGREE SEQUENCE D1 : ", conf->degree1, " D2 : ", conf->degree2);
    double N1 =  0.5 * conf->n_vertex; // Nodi di grado DEGREE 1
    double N2 =  0.5 * conf->n_vertex; // Nodi di grado DEGREE 2

    double P1 = 1.0 / N1;
    double P2 = 1.0 / N2;

    /* Probabilità che una pallina venga assegnata a un nodo di classe 1 */
    double Q1 = (0.5 * conf->degree1) / (0.5 * conf->degree1 + 0.5 * conf->degree2);

    /* Probabilità che una pallina venga assegnata a un nodo di classe 2*/
    double Q2 = (0.5 * conf->degree2) / (0.5 * conf->degree1 + 0.5 * conf->degree2);

    for(int j = 0; j <= conf->timer; j ++) vector_pdf[j] = 0;

    for(int k1 = 0; k1 <= conf->timer; k1++) {
        double t = gsl_ran_binomial_pdf(k1, Q1, conf->timer);
        for(int i = 0; i <= k1; i++) {
            vector_pdf[i] += t * (0.5 * gsl_ran_binomial_pdf(i, P1, k1));
        }
        for(int i = 0; i <= conf->timer - k1; i++) {
            vector_pdf[i] +=  t * (0.5 * gsl_ran_binomial_pdf(i, P2, conf->timer - k1));
        }
    }

    //CHECK IF SUM OF PDF EQUALS TO 1.0
    double somma = 0.0;
    for(int i = 0; i <= conf->timer; i++) {
        somma += vector_pdf[i];
        fprintf(output, "%d\t%f\n", i, vector_pdf[i]);
    }
    printf("Somma vector_pdf %f\n", somma);
    
}


/* Calcola e scrive su file la distribuzione uniforme teorica
    1° colonna : (asse X) n palline
    2° colonna : (asse Y) frequenza per il n palline

*/
void write_theoric(TYPE_OF_THEORIC type) {

    double *vector_pdf = (double *) malloc(sizeof(double) * (conf->timer + 1));
    double *vector_cdf = (double *) malloc(sizeof(double) * conf->timer + 1);

    get_file_name(THEORIC);
    // strcat(directory, file_name);
    char theoric_suffix[128];
    
    switch(type) {
        case BINOMIAL:
            /* Binomial distribution*/
            write_theoric_bin(vector_pdf, vector_cdf);
            strcpy(theoric_suffix, "bin_");
            break;
        case GEOMETRIC:
            /* Geometric distribution */
            write_theoric_geo(vector_pdf, vector_cdf);
            strcpy(theoric_suffix, "geo_");
        case GEOM_GEOM_1:
            /* Particular geometric distribution:
                (DeathPhase 2 no WITH_STEP)    
            */
            write_theoric_geom_geom_1(vector_pdf, vector_cdf);
            strcpy(theoric_suffix, "geom_geom_1");
            break;
        case WITH_CHOOSE:
            /* Particular theoric curve 
                (DeathPhase2 and WITH_STEP)
            */
            write_theoric_with_choose(vector_pdf, vector_cdf);
            strcpy(theoric_suffix, "with_choose_");
        break;  
        case STABLE_DEGREE:
            /* Theoric and frequency of only stable nodes */
            write_theoric_stable_nodes(vector_pdf, vector_cdf);
            strcpy(theoric_suffix, "stable_nodes_");

        break;
    }
    
    strcat(theoric_suffix, file_name);
    strcat(directory, theoric_suffix);
    strcpy(conf->file_name_theoric, directory);

    FILE *output = fopen(directory, "w");
    write_parameters(output);
    
    // Write on file
    for(int i = 0; i <= conf->timer; i++)
        fprintf(output, "%d\t%f\t%f\n", i, vector_pdf[i], vector_cdf[i]);
    printf("FILE NAME THEORIC : %s\n", directory);
}

