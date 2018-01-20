set terminal png enhanced  size 1200, 900
#set output 'radden.thymine.png'
#set output 'radden.thymine_polarized.png'
set output 'radden.thymine_5water_bigger.radius.png'  
#set key center top title " "
#set title ""
#set xlabel "Radius [a.u]"
#set ylabel " []"
#
set xtics 10,0.5,20
#set xrange [12.5: 15.5]
set yrange [0:0.0025]
#plot for [i=1:91] 'plots_thymine/density.matrix.'.i using 1:3 smooth unique with lp title 'density matrix '.i
plot for [i=1:91] 'plots_thymine_5water_bigger.radius/density.matrix.'.i using 1:3 smooth unique with lp title 'density matrix '.i
#
#do for [i=1:15] {
#  set terminal png enhanced  size 1200, 900
#  set output 'density.matrix.'.i.'.png'
#  set title ""
#  set xlabel ""              
#  set ylabel ""
#  unset key                               
#  plot 'density.matrix.'.i. using 1:3 with line pt 5
#}
set output   
set terminal x11
#
#plot 'density.matrix.1'u 1:3, 'density.matrix.2' u 1:3, 'density.matrix.3' using 1:3, 'density.matrix.4' using 1:3, 'density.matrix.5' using 1:3, 'density.matrix.6' using 1:3, 'density.matrix.7' using 1:3, 'density.matrix.8' using 1:3, 'density.matrix.9' using 1:3, 'density.matrix.10' using 1:3, 'density.matrix.11' using 1:3, 'density.matrix.12' using 1:3, 'density.matrix.13' using 1:3, 'density.matrix.14' using 1:3, 'density.matrix.15' using 1:3
