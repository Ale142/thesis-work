reset
set pointsize 1
set grid xtics
set grid ytics

set title "GRAFICO FREQUENZE SIMULAZIONE"
set xlabel "N PALLINE"

# Frequence medie e singola RUN
set term qt 0
set logscale y 
# set xrange [0:100]
plot "results/regular_graph/frequency/avg_regular_graph_frequency_1000000.dat" u 1:2 with boxes

# set term qt 1
# set xrange [0:1000]
# plot "avg_simulazione_10000.dat" u 1:2 with points,  "simulazione_10000.dat" u 1:2 with points,\
# "theoric_10000.dat" u 1:2 with points

# set term qt 2
#set xrange [0:10000]
# plot "avg_simulazione_100000.dat" u 1:2 with points,  "simulazione_100000.dat" u 1:2 with points,\
# "theoric_100000.dat" u 1:2 with points


