set title "Approximate Pressure."
set ylabel "p"
set xlabel "Length of Shock Tube (km)"
plot "ShockType1_Method1_Neutrinos1_ReconType1_Res1600_t00001375.csv" using 2:10 title "Approximate Solution"
set output "p.eps"
set terminal eps
replot