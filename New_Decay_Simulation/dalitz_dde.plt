set terminal qt
set title "Dalitz Plot for DDE Decay"
#set nokey
#set grid
#set xlabel "Excitation Energy"
#set ylabel "Efficiency"
set xrange[-0.5:0.5]
set yrange[-0.5:0.5]
m="dalitz_dde.txt"
plot m using 1:2 with points lc rgb "purple"
