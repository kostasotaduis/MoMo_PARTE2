pi = 4.0*atan(1.0)
mb(x) = 4 * pi * (2 * pi * 100.0)**(-1.5) * x*x * exp(-(x*x)/(2*100.0))



bw = 2.0
set boxwidth bw
b(x) = bw * floor(x/bw) + 0.5*bw
stats "vel0.dat" nooutput
N1 = STATS_records
stats "Velocidades.dat" nooutput
N2 = STATS_records
set xlabel "|v|"
set style fill solid
plot [0:30] "vel0.dat" u (b($1)):(1.0/(bw*N1*5)) smooth freq w boxes fc "#c63637" ti "Bimodal", \
"Velocidades.dat" u (b($1)):(1.0/(bw*N2)) smooth freq w boxes fc "#885083c1" ti "Equilibrium", mb(x) w l lw 3 lc "#5083c1" ti "Maxwell-Boltzmann"