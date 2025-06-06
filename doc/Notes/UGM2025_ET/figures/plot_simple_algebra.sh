#!/bin/sh
cd "${0%/*}" || exit                        # Run from this directory
#------------------------------------------------------------------------------

# settings

    # operand setups
    setups="
    DFSEM
    DFM
    FSM
    "


#------------------------------------------------------------------------------

image="simple_algebra.png"

gnuplot<<simple_algebra
set terminal pngcairo # font "helvetica,20" size 5000, 1000
set output "$image"
set xlabel "nCells"
set logscale x
set ylabel "Time [s]"
set logscale y
set grid
#set xrange [0:10000000]
#set yrange [0:1]
#set key left top
#set key samplen 2
#set key spacing 0.75
#set multiplot layout 1,11 title "Setup: $setup" noenhanced

plot \
"simple_algebra.txt" using 1:2 title "Intermediate field" w l, \
"simple_algebra.txt" using 1:3 title "Expression templates" w l

simple_algebra

#------------------------------------------------------------------------------
