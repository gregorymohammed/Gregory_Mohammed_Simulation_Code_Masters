set title "Approximate Velocity."
set ylabel "v"
set xlabel "Length of Shock Tube (km)"
plot "ShockType1_Method1_Neutrinos1_ReconType1_Res1600_t00001375.csv" using 2:7 title "Approximate Solution"
set output "v.eps"
set terminal eps
replot