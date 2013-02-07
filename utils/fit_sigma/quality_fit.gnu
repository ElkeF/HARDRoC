set output "/Users/elkefasshauer/Documents/quality_check/quality.ps"
set terminal postscript enhanced color solid "Helvetica" 20

# Titel, Axen Beschriftungen
set xrange [13:23]
set xlabel "R [angstrom]"
set ylabel "TRDM"

factor = 46.461521
const  = -5.8E-5

f(x) = factor * exp(-alpha*x) + const

plot "../../exppar_lib/sigma/Xe.dat" using 1:2 with points #, \
#     f(x) with lines
