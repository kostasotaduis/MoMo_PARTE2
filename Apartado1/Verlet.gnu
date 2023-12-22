set key center right
set xlabel "t (s)"
set ylabel "E (kJ/mol)"
plot [0:2] "p4_dt0000001.dat" u 1:2 w l lc "#5083c1" ti "E_{k}",\
 "" u 1:3 w l lc "#c63637" ti "E_{p}", \
"" u 1:($2+$3) w l lc "#000000" lw 2 ti "E_{tot}"

