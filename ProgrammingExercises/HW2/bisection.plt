reset
set xlabel 'N'
set ylabel 'iterations'
set log xy
set terminal pdf
set output 'asgn2_bisection.pdf'
plot 'bisection.dat' using 1:2 title 'sqrt(2)' with line, \
'bisection.dat' using 1:3 title 'cbrt(2)' with line, \
set terminal wxt
set output
exit