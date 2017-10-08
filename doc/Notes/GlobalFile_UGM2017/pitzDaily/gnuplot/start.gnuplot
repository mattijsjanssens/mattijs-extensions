load "settings.gnuplot"
set output "start.eps"
set title "Startup"
set yrange 0:400
plot \
    "uncollated_noTimeStampMaster/start.txt" title "Uncollated" \
    with lines lt 1, \
    "uncollated_timeStampMaster/start.txt" title "Uncollated+timeStampMaster" \
    with lines lt 2, \
    "collated_noTimeStampMaster/start.txt" title "Collated" \
    with lines lt 3, \
    "collated_timeStampMaster/start.txt" title "Collated+timeStampMaster" \
    with lines lt 4 \
