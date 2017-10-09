load "settings.gnuplot"
set output "write.eps"
set title "Writing"
set yrange [0:150]
plot \
    "uncollated_noTimeStampMaster/write.txt" title "Uncollated" \
    with linespoints lt 1, \
    "uncollated_timeStampMaster/write.txt" title "Uncollated+timeStampMaster" \
    with linespoints lt 2, \
    "collated_noTimeStampMaster/write.txt" title "Collated" \
    with linespoints lt 3, \
    "collated_timeStampMaster/write.txt" title "Collated+timeStampMaster" \
    with linespoints lt 4 \
