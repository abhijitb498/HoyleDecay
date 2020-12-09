# theta_phi_DDL.plt
set terminal x11 0
set title "Decay in Linear Chain"
#set nokey
#set grid
set xlabel "theta"
set ylabel "phi"
m="theta_phi_DDL.txt"
plot m using 1:2 with points lc rgb "green"
