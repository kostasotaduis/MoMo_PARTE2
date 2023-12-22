set key top right
set xlabel "t (s)"
set ylabel "p (kg*m/s)"
plot [0:2][100:106] "p4_dt000001.dat" u 1:4 w l lc "#5083c1" lw 2 noti

