reset
set pointsize 1
set grid xtics
set grid ytics

set title "# WATTS-STROGATZ small-world, n_vertex: 100, N_RAN: 100, TIMER: 100000, DEATH_PROB: 0.015000, DEATH_PHASE: 2, K: 1 REWIRE_P = 0.0, DIM = 1, NEI = 10 LOOPS = True MULTIPLE = True"
set xlabel "N PALLINE"
set style line 3 lw 3 lc rgb '#0044a5' ps 2 pt 9 pi 1

# Frequence medie e singola RUN
set term qt 0 
set xrange[:50]
# plot "results/regular_graph/frequency/avg_regular_graph_frequency_1000000.dat" u 1:3 with points t "Simulazione (frequenze medie con n_ran)", "results/regular_graph/theoric/geom_geom_1regular_graph_theoric_1000000.dat" u 1:2 with points t "GEOM/GEOM/1"
# plot "results/regular_graph/frequency/step_avg_regular_graph_frequency_1000000.dat" u 1:3 with boxes t "Frequency"
plot "results/watts_strogatz/frequency/step_avg_watts_strogats_frequency_100000.dat" u 1:3 with boxes t "Frequency", "results/watts_strogatz/theoric/with_choose_watts_strogatz_100000.dat" u 1:2 with boxes t "Theoric"


# set term qt 1
# plot "results/regular_graph/theoric/with_choose_regular_graph_theoric_1000000.dat" u 1:2 with boxes t "Theoric"

# set term qt 1
# set xrange[:50]
# plot "results/regular_graph/theoric/with_choose_regular_graph_theoric_10000.dat" u 1:2 with boxes


# set term qt 1
# set yrange [:0.0005]
# plot "avg_simulazione_10000.dat" u 1:3 with boxes,  "simulazione_10000.dat" u 1:3 with p

# set term qt 2
# set yrange [:0.000040]
# plot "avg_simulazione_100000.dat" u 1:3 with lp,  "simulazione_100000.dat" u 1:3 with p