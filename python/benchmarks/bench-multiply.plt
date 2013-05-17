a = 1; b = 2.8; c = 1; d = 3
FFPACK_f(x) = a*x**b
NZMATH_f(x) = c*x**d
fit FFPACK_f(x) "./multiply.dat" using 1:2 via a,b
fit NZMATH_f(x) "./multiply.dat" using 1:4 via c,d
set logscale y
set terminal postscript eps color font "" 18
set output "multiply.eps"
set key bottom right
set xlabel "Size of matrices (n-by-n)"
set ylabel "Running time (s)"
plot "multiply.dat" using 1:2:3 lt 1 title "FFPACK" with errorbars, FFPACK_f(x) title sprintf("%.3g*x**%.3g", a, b) lt 1, "multiply.dat" using 1:4:5 lt 2 title "NZMATH" with errorbars, NZMATH_f(x) title sprintf("%.3g*x**%.3g", c, d) lt 2
