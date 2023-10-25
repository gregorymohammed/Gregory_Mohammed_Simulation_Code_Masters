set title "Approximate Rest Density."
set ylabel "rho0"
set xlabel "Length of Shock Tube (km)"
plot "ShockType1_Method1_Neutrinos1_ReconType1_Res1600_t00001375.csv" using 2:8 title "Approximate Solution"
set output "rho0.eps"
set terminal eps
replot