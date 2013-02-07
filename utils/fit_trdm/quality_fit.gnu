set output "/Users/elkefasshauer/Documents/quality_check/quality.ps"
set terminal postscript enhanced color solid "Helvetica" 20

# Titel, Axen Beschriftungen
set xrange [3.5:14]
set xlabel "R [angstrom]"
set ylabel "TRDM"

factor = 46.461521
alpha  = 1.734922
const  = -5.8E-5

f(x) = factor * exp(-alpha*x) + const

plot "ArXe.txt" using 1:6 with points, \
     f(x) with lines
