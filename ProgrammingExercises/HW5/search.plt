reset
set xlabel 'N'
set ylabel 'seconds'
set log xy
set terminal pdf
set output 'search.pdf'
plot 'search.dat' using 1:2 title 'linear' with line, \
'search.dat' using 1:3 title 'binary' with line, \
'search.dat' using 1:4 title 'dictonary' with line
exit