epsilon = 0.998
tE = epsilon
sigma = 3.4e-10
m = 40 * (6.022e23)**(-1)

tp = (epsilon * m)/(sigma**3) * (1.0/1.0e6)

tt = ((m*1e-3)/(epsilon*1e3))**(0.5) * sigma * 1e-10
trho = (m)/((sigma)**3) * (1.0/1.0e6)
set key center right
set ylabel "p (Pa)"
set xlabel "{/Symbol r} (g/cm^{3})"
plot "Evsrho.txt" u ($2*trho):($5*tp) w lp pt 7 lc "#42ab49" noti