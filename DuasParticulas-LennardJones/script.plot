#/bin/gnuplot
set title "Energias vs tempo"
set xrange[0:5000]
set yrange[0:7]
set ylabel "Energia"
set xlabel "T"
 
plot "energia.dat"  u 1:2 w lp title "energia cinetica",
replot "energia.dat" title "energia potencial", u 1:3 w lp
replot "energia.dat" title "energia total", u 1:4 w lp

set terminal jpeg
set output "energias.jpeg"
replot
