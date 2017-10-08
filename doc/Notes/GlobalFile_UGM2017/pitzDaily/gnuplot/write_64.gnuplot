load "settings.gnuplot"
set output "write_64.eps"
set title "Writing"
set xtics 10
set yrange [0:40]
plot \
    "collated_timeStampMaster/write_64.txt" title "Collated+timeStampMaster" \
    with linespoints lt 1 \
