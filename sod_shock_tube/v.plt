set ylabel "The Fluid Velocity"
set xlabel "Length of Shock Tube (km)"
plot "ShockType1_Method1_Neutrinos0_ReconType1_Res1600_t00002150.csv" using 2:7 title "Approximate Solution", "Riemann_Exact_Solution.csv" using 1:2 with lines title "Exact Solution", "Riemann_Exact_Solution.csv" using 1:8 with linespoints title "Difference"
set output "v.eps"
set terminal eps
replot