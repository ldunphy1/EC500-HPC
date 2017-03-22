reset
set xlabel 'N'
set ylabel 'seconds'
set log xy
set terminal pdf
set output 'jacobi_timing.pdf'
plot 'my_jacobi_mpi_timing_4.dat' using 1:2 title '4 processes' with line, \
'my_jacobi_mpi_timing_8.dat' using 1:2 title '8 processes' with line, \
'my_jacobi_mpi_timing_16.dat' using 1:2 title '16 processes' with line
exit