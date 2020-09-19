#ifndef __FREQUENCY__
#define __FREQUENCY__
#include <stdio.h>
#include <unistd.h>
#include <string.h>

void count_frequency(int *frequency, int *palline);
void count_k_balls(int time);
void print_frequency(int *frequency);
void write_frequency(int *frequency);
void write_k_balls();

void calculate_average_frequency(double *average_frequency, double **frequency_nran);
void write_average_frequency(double *avg_frequency);
void write_single_frequency(double *avg_freq_degree1, double *avg_freq_degree2);
void count_frequency_single_degree(double *frequency, int *palline, int DEGREE);
#endif