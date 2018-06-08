
set yrange [-0.01:0.09]
set xrange [-0.1:2.1]
set terminal png size 1366,768
set output "A_Wykres3_T.png"
set title 'Wykres zalezność maksymalnej wartosci bezwzglednej bledu w funkcji czasu t'
set ylabel 'y'
set xlabel 'x'

set grid 
plot 'Wykres3_laasonen_thomasa.txt'  using 1:2 with lines title "Thomasa",\
 #'Wykres3_laasonen_SOR.txt'  using 1:2  with lines title  "Numeryczne",\
 


set yrange [-0.01:0.09]
set xrange [-0.1:2.1]
set terminal png size 1366,768
set output "A_Wykres3_S.png"
set title 'Wykres zalezność maksymalnej wartosci bezwzglednej bledu w funkcji czasu t'
set ylabel 'y'
set xlabel 'x'

set grid 
plot 'Wykres3_laasonen_SOR.txt'  using 1:2  with lines title  "SOR",\
 

