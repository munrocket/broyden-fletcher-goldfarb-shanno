#dfp isolines
set table "penalty_iso.dat"
set contour
unset surface
set view 0, 0
set xrange [-2.1:-.4]
set yrange [-1:1]
set isosample 50
set cntrparam bspline
set cntrparam levels auto 60
splot x**2 + y**2

#plot
set term postscript eps
set output "penalty.eps"
#set xtics 0.2
#set ytics 0.2
#set mxtics 2
#set mytics 2
unset table
unset key
plot 	"penalty_iso.dat" using 1:2 with line linetype 0, \
	"points1.dat" u 1:2 w linespoints lt 1 pointtype

#plot in eps
set term window
replot
pause -1