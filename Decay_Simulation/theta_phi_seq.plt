# theta_phi_seq.plt
#set terminal png
set terminal x11 0
#set output "theta_phi_seq.png"
set title "Sequential Decay"
#set nokey
#set grid
set xlabel "theta"
set ylabel "phi"
m="theta_phi_seq.txt"
plot m using 1:2 with points lc rgb "red"
