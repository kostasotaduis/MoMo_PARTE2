set key top left
set xlabel "t (s)"
set ylabel "E (kJ/mol)"
plot [0:2][-5000:80000] "p4_Eulerdt00001.dat" u 1:2 w l lc "#5083c1" ti "E_{k}",\
 "" u 1:3 w l lc "#c63637" ti "E_{p}", \
"" u 1:($2+$3) w l lc "#000000" lw 2 ti "E_{tot}"

