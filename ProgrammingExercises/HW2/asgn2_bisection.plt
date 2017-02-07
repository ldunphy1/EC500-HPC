reset
f(x)=c0 + c1*log(x) + c2*log(log(x))
c0 = 1; x1 = 1; x2 = 1;
fit f(x) 'bisection.dat' using 1:2 via c0, c1, c2
set xlabel 'N'
set ylabel 'iterations'
set log xy
set key top left
set terminal pdf
set output 'asgn2_bisection.pdf'
plot 'bisection.dat' using 1:2 title 'sqrt(2)' with line, \
'bisection.dat' using 1:3 title 'cbrt(2)' with line, \
c0 + c1*log(x) + c2*log(log(x)) title 'fitted curve'
set terminal wxt
set output
exit