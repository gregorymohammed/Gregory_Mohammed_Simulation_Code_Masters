set title "The Initial Energy."
set ylabel "energy"
set xlabel "Length of Shock Tube (km)"
plot "ShockType1_Method1_Neutrinos0_ReconType1_Res1600_t00000000.csv" using 2:9 title "Energy"
set output "initialEps.eps"
set terminal eps
replot