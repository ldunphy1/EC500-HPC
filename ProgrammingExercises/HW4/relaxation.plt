reset
set xlabel 'N'
set ylabel 'iterations'
set terminal pdf
set output 'relaxation.pdf'
plot 'relaxation.dat' using 1:2 title 'jacobi(N)' with line
exit