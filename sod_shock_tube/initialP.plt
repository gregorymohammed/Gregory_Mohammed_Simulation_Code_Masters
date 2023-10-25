set title "The Initial Pressure."
set ylabel "pressure"
set xlabel "Length of Shock Tube (km)"
plot "ShockType1_Method1_Neutrinos0_ReconType1_Res1600_t00000000.csv" using 2:10 title "Pressure"
set output "initialP.eps"
set terminal eps
replot