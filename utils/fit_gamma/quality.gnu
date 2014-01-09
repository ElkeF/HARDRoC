#set output "/Users/elkefasshauer/Documents/quality_check/quality.ps"
#set terminal postscript enhanced color solid "Helvetica" 20

# Titel, Axen Beschriftungen
set logscale xy
set xrange [2:23]
set yrange [1.0e-7:10]
set xlabel "R [AA]"
set ylabel "Gamma [eV]"

#Ne2
alpha   =  19.643762
beta    =   2.789916
gamma   =   6.204031 # for R< 4.00 r^-6 term

#NeAr
#gamma   =            23.797274 # for fit only distances above eq are used
#alpha   = 1211057096516.005371 # afterwards add alpha, beta fit
#beta    =            11.074016

f(x) = alpha * exp(-beta*x) + gamma/x**6
#f(x) = gamma/x**6


plot "Ne2.dat" using 1:2 with points , \
     f(x) with lines
