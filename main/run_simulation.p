reset
set pointsize 1
set grid xtics
set grid ytics
set title'# REGULAR GRAPH, conf->n_vertex : 100, N_RAN: 1, TIMER: 10000,T_TRANS: 0, conf->degree1: 10, DEATH_PROB: 0.012500, DEATH_PHASE: 2 K: 3'
set xrange[:50]
plot 'results/regular_graph/frequency/step_avg_regular_graph_frequency_10000.dat' u 1:3 with boxes t 'Frequency', 'results/regular_graph/theoric/with_choose_regular_graph_theoric_10000.dat' u 1:2 with boxes t 'Theoric'
