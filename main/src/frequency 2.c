#include "../lib/frequency.h"
#include "../lib/simulation.h"


// void count_frequency_transitory(int *frequency, int *palline) {

// }

void count_frequency(int *frequency, int *palline) {
    int count = 0;
    for(int i = 0; i <= conf->timer; i++) {

        for(int j = 0; j < conf->n_vertex; j++) {
            // In quante urne trovo i palline?
            if(i == palline[j])
                count++;
        }
        // printf("Count of %d balls: %d\n", i, count);
        frequency[i] += count;
        count = 0;
    }
    // printf("\n");
}

void count_k_balls(int time) {
    int count = 0;
    for(int i = 0; i < conf->n_vertex; i++)
        if(palline[i] >= 2)
           count++;
    zero_balls[time] = count + 1;
}

void write_k_balls() {
    FILE *output = fopen("zero_balls.dat", "w");
    fprintf(output, "# URNE CON 0 PALLINE IN FUNZIONE DEL TEMPO");
    for(int i = 0; i <= conf->timer; i++) 
        fprintf(output, "%d\t%d\n", i, zero_balls[i]);
    fclose(output);
}

void print_frequency(int *frequency) {
    for(int i = 0; i < conf->timer; i++)
        printf("(%d - %d) \t", i, frequency[i]);
}

/* Scrivo su file le frequenze calcolate:
    1° colonna : (asse X) n palline
    2° colonna : (asse Y) frequenza calcolata per il n palline
*/
void write_frequency(int *frequency) {
    get_file_name(FREQUENCY);

    strcat(directory, file_name);
    strcpy(conf->file_name_frequency, directory);
    printf("File name %s\n", directory);

    // strcat(directory, file_name);
    FILE *output_file = fopen(directory, "w");
    
    write_parameters(output_file);
    fprintf(output_file, "%s\t%s\t%s\n","# N_palline(X)", "frequency(Y1)", "probability(Y2)");
    for(int i = 0;i <= conf->timer; i++){
        if(frequency[i] != 0)
            fprintf(output_file, "%d\t%d\t%f\n", i, frequency[i], (double)frequency[i] / conf->timer);

    }
}

/* Calcola la frequenza media */
void calculate_average_frequency(double *average_frequency, double **frequency_nran) {
    double somma = 0.0;
    for(int i = 0;i <= conf->timer; i++) {
        for(int j = 0; j < conf->n_ran; j++) {
            somma += frequency_nran[j][i];
        }
        average_frequency[i] = somma / conf->n_ran;
        somma = 0.0;
    }
}

void write_average_frequency(double *avg_frequency) {
    get_file_name(FREQUENCY);

    strcat(directory, file_name);
    strcpy(conf->file_name_frequency, directory);

    printf("FILE NAME AVG FREQUENCY : %s\n", directory);
    // printf("file name %s\n", directory);
    FILE *output_file = fopen(directory, "w");
    write_parameters(output_file);
    fprintf(output_file, "%s\t%s\t%s\n","# N_palline(X)", "AVG frequency(Y1)", "probability(Y2)");
    
    double *check_pdf = (double *) malloc(sizeof(double) * (conf->timer + 1));
    for(int i = 0;i <= conf->timer; i++){
        check_pdf[i] = 0;
        if(avg_frequency[i] != 0) {
            check_pdf[i] = (avg_frequency[i]) /  conf->n_vertex;
            fprintf(output_file, "%d\t%f\t%f\n", i, avg_frequency[i], ((double)avg_frequency[i]) / conf->n_vertex);
        }
    }

    double somma = 0;
    for(int i = 0; i <= conf->timer; i++) {
        somma += check_pdf[i];
    }
    printf("SUM_PDF_FREQUENCY %f\n", somma);
    fclose(output_file);
}


/* Usata nella simulazione grafo Degree_sequence e due gradi Degee1, Degree2.
   Scrive su file le frequenze medie per i due gradi
 */
void write_single_frequency(double *avg_freq_degree1, double *avg_freq_degree2) {
    get_file_name(FREQUENCY);

    double N1 =  0.5 * conf->n_vertex; // Nodi di grado DEGREE 1
    double N2 =  0.5 * conf->n_vertex; // Nodi di grado DEGREE 2

    char *suffix, *path_degree1, *path_degree2;
    suffix = (char *)malloc(sizeof(char) * 100);
    path_degree1 = (char *) malloc(sizeof(char) * 100);
    path_degree2 = (char *) malloc(sizeof(char) * 100);
    
    strcpy(path_degree1, directory);
    strcpy(path_degree2, directory);

    strcpy(suffix, "degree1_");
    strcat(suffix, file_name);
    strcat(path_degree1, suffix);

    strcpy(suffix, "degree2_");
    strcat(suffix, file_name);
    strcat(path_degree2, suffix);

    FILE *file_degree1 = fopen(path_degree1, "w");
    FILE *file_degree2 = fopen(path_degree2, "w");

    write_parameters(file_degree1);
    write_parameters(file_degree2);

    fprintf(file_degree1, "%s\t%s\t%s\n","# N_palline(X)", "AVG frequency_degree1(Y1)", "probability_degree1(Y2)");
    fprintf(file_degree2, "%s\t%s\t%s\n","# N_palline(X)", "AVG frequency_degree2(Y1)", "probability_degree2(Y2)");


    for(int i = 0;i <= conf->timer; i++){
      
        if(avg_freq_degree1[i] != 0) {
            fprintf(file_degree1, "%d\t%f\t%f\n", i, avg_freq_degree1[i], avg_freq_degree1[i] / N1);
        }
        if(avg_freq_degree2[i] != 0) {
            fprintf(file_degree2, "%d\t%f\t%f\n", i, avg_freq_degree2[i], avg_freq_degree2[i] / N2 );
        }
    }

    printf("WROTE SINGLE FREQUENCY IN FILE : %s\n", directory);
    fclose(file_degree1);
    fclose(file_degree2);
}

/* Utilizzata per calcolare le frequenze medie per un particolare grado */
void count_frequency_single_degree(double *frequency_single, int *palline, int DEGREE) {
    int count = 0;
    for(int i = 0; i <= conf->timer; i++) {

        for(int j = 0; j < conf->n_vertex; j++) {
            // In quante urne trovo i palline?
            
            if(get_degree(j) == DEGREE) {
                // printf("Vertex %d, Degree %d\n", j, DEGREE);
                if(i == palline[j])
                    count++;
            }
        }
        frequency_single[i] = count;
        count = 0;
    }
}