set output "/Users/elkefasshauer/Documents/quality_check/quality.ps"
set terminal postscript enhanced color solid "Helvetica" 20

# Titel, Axen Beschriftungen
set xrange [13:23]
set xlabel "Energy [eV]"
set ylabel "TRDM"

quad    =     0.695058
lin     =   -41.741611
const   =   773.748197
oneover = -3655.360124

f(x) = quad*x**2 + lin*x + const + oneover/x

plot "../../exppar_lib/sigma/Xe.dat" using 1:2 with points , \
     f(x) with lines
