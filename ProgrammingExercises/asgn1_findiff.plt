reset
set xlabel 'x'
set ylabel 'relative error'
set log x
set log y
set format x '%.2e'; set xtics font ',7'
set format y '%.2e'; set ytics font ',7'
set terminal pdf
set output 'asgn1_findiff.pdf'
plot 'asgn1_findiff_cppOut.txt' using 1:(abs(($2-$5)/$5)) title 'forward' with line, \
'asgn1_findiff_cppOut.txt' using 1:(abs(($3-$5)/$5)) title 'backward' with line, \
'asgn1_findiff_cppOut.txt' using 1:(abs(($4-$5)/$5)) title 'central' with line
set terminal wxt
set output
exit