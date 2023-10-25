set ylabel "u"
set xlabel "Length of Shock Tube (km)"
plot "ShockType1_Method1_Neutrinos0_ReconType1_Res800_t00001120.csv" using 2:9 title "Approximate Solution", "Riemann_Exact_Solution.csv" using 1:4 with lines title "Exact Solution", "Riemann_Exact_Solution.csv" using 1:9 with linespoints title "Difference"
set output "u.eps"
set terminal eps
replot