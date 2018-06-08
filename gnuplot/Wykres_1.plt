set yrange [-4.5:-1]
set xrange [-1.7:0.2]
#set terminal wxt size 1280,900
set terminal png size 1366,768
set output "A_W1_thomas.png"
set title 'Wykres zaleznosci maksymalnejwartosci bezwzglednaj bledu obserwowanej dla t_max, w funkcji kroku przestrzennegp h (w skali logarytmicznej) '


set ylabel 'log10(|bledu|)'
set xlabel 'log10(|kroku|)
set grid 

plot "W1_bmax_laasonen_thomasa_log.csv" using 1:2 with lines title "Thomas",\



set yrange [-4.5:-1]
set xrange [-1.7:0.2]
#set terminal wxt size 1280,900
set terminal png size 1366,768
set output "A_W1_SOR.png"
set title 'Wykres zaleznosci maksymalnejwartosci bezwzglednaj bledu obserwowanej dla t_max, w funkcji kroku przestrzennegp h (w skali logarytmicznej) '


set ylabel 'log10(|bledu|)'
set xlabel 'log10(|kroku|)
set grid 

plot "W1_bmax_laasonen_SOR_log.csv" using 1:2 with lines title "SOR",\




