set logscale y
set xrange [0:51]
set terminal postscript eps color font "" 18
set output "generator.eps"
set key bottom right
set xlabel "Dimension of algebra"
set ylabel "Running time (s)"
plot "generator-mean-stdev.dat" using 1:2:3 lt 1 with errorbars title "GF(2)", "generator-7-mean-stdev.dat" using 1:2:3 lt 2 with errorbars title "GF(7)"
