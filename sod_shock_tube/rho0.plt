set ylabel "The Rest Density (rho0) in Sod Units"
set xlabel "Length of Shock Tube (km)"
plot "ShockType1_Method1_Neutrinos0_ReconType1_Res1600_t00002150.csv" using 2:8 title "Approximate Solution", "Riemann_Exact_Solution.csv" using 1:3 with lines title "Exact Solution", "Riemann_Exact_Solution.csv" using 1:6 with linespoints title "Difference"
set output "rho0.eps"
set terminal eps
replot