# theta_phi_DDphi.plt
set terminal x11 0
set title "Decay in Equal Phase Space"
#set nokey
#set grid
set xlabel "theta"
set ylabel "phi"
m="theta_phi_DDphi.txt"
plot m using 1:2 with points lc rgb "magenta"
