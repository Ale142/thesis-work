reset
set pointsize 1
set grid xtics
set grid ytics

set title "GRAFICO THEORIC"
set xlabel "N PALLINE"


# Frequenze medie PDF e confronto con curva teorica
reset
# set term qt 1

# set term qt 0
set logscale y
# set yrange [:0.15]
plot "results/regular_graph/frequency/avg_regular_graph_frequency_100000.dat" u 1:3 with boxes, "results/regular_graph/theoric/geo_regular_graph_theoric_100000.dat" u 1:2 with boxes

# reset
# set term qt 1
# set xrange [0:200]
# plot "avg_simulazione_10000.dat" u 1:3 with points, "theoric_10000.dat" u 1:2 with points

# reset 
# set xrange [0:1500]
# set term qt 2
# plot "avg_simulazione_100000.dat" u 1:3 with points, "theoric_100000.dat" u 1:2 with points