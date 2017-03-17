reset
set xlabel 'N'
set ylabel 'residual'
set terminal pdf
set output 'residual.pdf'
plot 'residual.dat' using 1:2 title 'residual' with line
exit