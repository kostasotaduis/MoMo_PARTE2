epsilon = 0.998
tE = epsilon
sigma = 3.4e-10
m = 40 * (6.022e23)**(-1)

tt = ((m*1e-3)/(epsilon*1e3))**(0.5) * sigma * 1e-10
trho = (m)/((sigma)**3) * (1.0/1.0e6)
set key center right
set ylabel "E (kJ/mol)"
set xlabel "{/Symbol r} (g/cm^{3})"
plot "Evsrho.txt" u ($2*trho):($4*tE) w lp pt 7 lc "#5083c1" ti "E_{k}", \
"" u ($2*trho):($3*tE) w lp pt 7 lc "#c63637" ti "E_{p}", \
"" u ($2*trho):(($3+$4)*tE) w lp pt 7 lc "#00000" ti "E_{tot}"