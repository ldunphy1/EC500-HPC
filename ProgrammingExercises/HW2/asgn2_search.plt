reset
set xlabel 'N'
set ylabel 'iterations'
set log xy
set key top left
set terminal pdf
set output 'asgn2_search.pdf'
plot 'search_timing.dat' using 1:2 title 'linear search' with line, \
'search_timing.dat' using 1:3 title 'binary search' with line, \
'search_timing.dat' using 1:4 title 'dictionary search' with line
set terminal wxt
set output
exit