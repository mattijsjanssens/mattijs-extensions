load "settings.gnuplot"
set output "start.eps"
set title "Startup"
plot \
    "uncollated_noTimeStampMaster/start.txt" title "Uncollated" \
    with linespoints lt 1, \
    "uncollated_timeStampMaster/start.txt" title "Uncollated+timeStampMaster" \
    with linespoints lt 2, \
    "collated_noTimeStampMaster/start.txt" title "Collated" \
    with linespoints lt 3, \
    "collated_timeStampMaster/start.txt" title "Collated+timeStampMaster" \
    with linespoints lt 4 \
