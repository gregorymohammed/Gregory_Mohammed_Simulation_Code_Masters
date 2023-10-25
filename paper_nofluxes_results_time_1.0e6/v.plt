set ylabel "v"
set xlabel "Length of Shock Tube (km)"
set key box
plot "ShockType1_Method1_Neutrinos1_ReconType1_Res1600_t00001734.csv" using 2:7 title "Time = 1.0e6 km"
set output "v.eps"
set terminal eps
replot