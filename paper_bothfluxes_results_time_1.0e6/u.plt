set ylabel "u"
set xlabel "Length of Shock Tube (km)"
set key box
plot "ShockType1_Method1_Neutrinos1_ReconType1_Res1600_t00001466.csv" using 2:9 title "Time = 1.0e6 km"
set output "u.eps"
set terminal eps
replot