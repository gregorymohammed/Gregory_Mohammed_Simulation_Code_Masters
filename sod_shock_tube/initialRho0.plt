set title "The Initial Rho0."
set ylabel "rho0"
set xlabel "Length of Shock Tube (km)"
plot "ShockType1_Method1_Neutrinos0_ReconType1_Res1600_t00000000.csv" using 2:8 title "Rest Mass"
set output "initialRho0.eps"
set terminal eps
replot