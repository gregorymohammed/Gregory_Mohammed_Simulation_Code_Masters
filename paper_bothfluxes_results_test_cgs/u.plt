set title "Approximate Internal Energy."
set ylabel "u"
set xlabel "Length of Shock Tube (km)"
plot "ShockType1_Method1_Neutrinos1_ReconType1_Res1600_t00001375.csv" using 2:9 title "Approximate Solution"
set output "u.eps"
set terminal eps
replot