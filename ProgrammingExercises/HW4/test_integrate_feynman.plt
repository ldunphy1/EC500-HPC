reset
set xlabel 'm'
set ylabel 'P(m)'
set terminal pdf
set output 'integral_scaling.pdf'
plot 'test_integrate_feynman.dat' using 1:2 title 'P(m)' with line, \
log(x) with line
exit