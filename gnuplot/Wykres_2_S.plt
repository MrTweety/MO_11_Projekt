
set yrange [-0.1:1.1]
set xrange [-10:10]
set terminal png size 1366,768
set output "A_SOR_Wykres2_i50.png"
set title 'Wykres rozwiązania analitycznego oraz numerycznego dla chwili t = 0.276726'
set ylabel 'y'
set xlabel 'x'

set grid 
plot 'rozwA50.txt'  using 2:3 with lines title "Analityczne",\
 'rozwSOR50.txt'  using 2:3 with points pt 8 lc 13 title "SOR",\
 

set yrange [-0.1:1.1]
set xrange [-10:10]
set terminal png size 1366,768
set output "A_SOR_Wykres2_i100.png"
set title 'Wykres rozwiązania analitycznego oraz numerycznego dla chwili t = 0.553452'
set ylabel 'y'
set xlabel 'x'

set grid 
set grid 
plot 'rozwA100.txt'  using 2:3 with lines title "Analityczne",\
 'rozwSOR100.txt'  using 2:3 with points pt 8 lc 13title "SOR",\
 
 
set yrange [-0.1:1.1]
set xrange [-10:10]
set terminal png size 1366,768
set output "A_SOR_Wykres2_i150.png"
set title 'Wykres rozwiązania analitycznego oraz numerycznego dla chwili t = 0.830178'
set ylabel 'y'
set xlabel 'x'

set grid 
plot 'rozwA150.txt'  using 2:3 with lines title "Analityczne",\
 'rozwSOR150.txt'  using 2:3 with points pt 8 lc 13 title "SOR",\
 

set yrange [-0.1:1.1]
set xrange [-10:10]
set terminal png size 1366,768
set output "A_SOR_Wykres2_i200.png"
set title 'Wykres rozwiązania analitycznego oraz numerycznego dla chwili t = 1.1069'
set ylabel 'y'
set xlabel 'x'

set grid 
plot 'rozwA200.txt'  using 2:3 with lines title "Analityczne",\
 'rozwSOR200.txt'  using 2:3 with points  pt 8 lc 13 title "SOR",\
 

set yrange [-0.1:1.1]
set xrange [-10:10]
set terminal png size 1366,768
set output "A_SOR_Wykres2_i250.png"
set title 'Wykres rozwiązania analitycznego oraz numerycznego dla chwili t = 1.38363'
set ylabel 'y'
set xlabel 'x'

set grid 
plot 'rozwA250.txt'  using 2:3 with lines title "Analityczne",\
 'rozwSOR250.txt'  using 2:3 with points  pt 8 lc 13 title "SOR",\
 

set yrange [-0.1:1.1]
set xrange [-10:10]
set terminal png size 1366,768
set output "A_SOR_Wykres2_i300.png"
set title 'Wykres rozwiązania analitycznego oraz numerycznego  dla chwili t = 1.66036'
set ylabel 'y'
set xlabel 'x'

set grid 
plot 'rozwA300.txt'  using 2:3 with lines title "Analityczne",\
 'rozwSOR300.txt'  using 2:3 with points  pt 8 lc 13 title "SOR",\
 

set yrange [-0.1:1.1]
set xrange [-10:10]
set terminal png size 1366,768
set output "A_SOR_Wykres2_i361.png"
set title 'Wykres rozwiązania analitycznego oraz numerycznego dla chwili t = 2.0'
set ylabel 'y'
set xlabel 'x'

set grid 
plot 'rozwA361.txt'  using 2:3 with lines title "Analityczne",\
 'rozwSOR361.txt'  using 2:3 with points pt 8 lc 13 title "SOR",\
 

