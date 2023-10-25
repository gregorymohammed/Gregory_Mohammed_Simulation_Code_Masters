set ylabel "The Pressure in Sod Units"
set xlabel "Length of Shock Tube (km)"
plot "ShockType1_Method1_Neutrinos0_ReconType1_Res1600_t00002150.csv" using 2:10 title "Approximate Solution", "Riemann_Exact_Solution.csv" using 1:5 with lines title "Exact Solution", "Riemann_Exact_Solution.csv" using 1:7 with linespoints title "Difference"
set output "p.eps"
set terminal eps
replot