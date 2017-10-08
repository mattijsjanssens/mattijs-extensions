load "settings.gnuplot"
set output "run.eps"
set title "Run 20 iterations"
set yrange [0:400]
plot \
    "uncollated_noTimeStampMaster/run.txt" title "Uncollated" \
    with linespoints lt 1, \
    "uncollated_timeStampMaster/run.txt" title "Uncollated+timeStampMaster" \
    with linespoints lt 2, \
    "collated_noTimeStampMaster/run.txt" title "Collated" \
    with linespoints lt 3, \
    "collated_timeStampMaster/run.txt" title "Collated+timeStampMaster" \
    with linespoints lt 4 \
