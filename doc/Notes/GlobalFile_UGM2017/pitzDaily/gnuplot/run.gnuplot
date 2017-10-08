load "settings.gnuplot"
set output "run.eps"
set title "Run 20 iterations"
plot \
    "uncollated_noTimeStampMaster/run.txt" title "Uncollated" \
    with lines lt 1, \
    "uncollated_timeStampMaster/run.txt" title "Uncollated+timeStampMaster" \
    with lines lt 2, \
    "collated_noTimeStampMaster/run.txt" title "Collated" \
    with lines lt 3, \
    "collated_timeStampMaster/run.txt" title "Collated+timeStampMaster" \
    with lines lt 4 \
