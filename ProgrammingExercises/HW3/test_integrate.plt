reset
set xlabel 'N'
set ylabel 'relative error'
set log y
set terminal pdf
set output 'x8Gauserr.pdf'
plot 'test_integrate.dat' using 1:(abs($5)+1e-20) title 'f(x)=x^8' with line
exit