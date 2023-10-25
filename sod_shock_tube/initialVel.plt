set title "The Initial Velocity."
set ylabel "velocity"
set xlabel "Length of Shock Tube (km)"
plot "ShockType1_Method1_Neutrinos0_ReconType1_Res1600_t00000000.csv" using 2:7 title "Velocity"
set output "initialVel.eps"
set terminal eps
replot