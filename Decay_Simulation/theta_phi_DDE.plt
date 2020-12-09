# theta_phi_DDE.plt
set terminal x11 0
set title "Decay with Equal Energy"
#set nokey
#set grid
set xlabel "theta"
set ylabel "phi"
m="theta_phi_DDE.txt"
plot m using 1:2 with points lc rgb "orange"
