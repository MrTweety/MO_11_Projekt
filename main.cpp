#include <iostream>
#include <cmath>
#include "projekt.h"
#include <unistd.h>
#include <ctime>

void czekaj( int iSekundy )
{
    for( clock_t koniec = clock() + iSekundy * CLOCKS_PER_SEC; clock() < koniec; )
        continue;

}


int main() {
    std::cout << "Hello, World!" << std::endl;

    int i=5;
    int ile_iteracji=pow(2.0,i+4);


    Projekt proj(ile_iteracji,0);
    proj.rozwiaz_analityczne();
    proj.rozwiaz_laasonen_thomasa();
    //proj.rozwiaz_laasonen_SOR();
    proj.save_gnuplot3(proj.rozwiazanieT,"a","a");
    proj.f_blad();
    proj.saveA(proj.blad_T,"bladT.csv");


//    Projekt proj(ile_iteracji,0);
//    proj.rozwiaz_analityczne();
//
//    time_t czasStart = time( NULL );
//    proj.rozwiaz_laasonen_thomasa();
//    time_t czasStop = time( NULL );
//    double czasT = difftime( czasStop, czasStart );
//    printf( "Uplynelo %.2fsek.\n", difftime( czasStop, czasStart ) );
//
//
//    czasStart = time( NULL );
//    proj.rozwiaz_laasonen_SOR();
//    czasStop = time( NULL );
//    double czasSOR = difftime( czasStop, czasStart );
//    printf( "Uplynelo %.2fsek.\n", difftime( czasStop, czasStart ) );
//
//    printf( "czasSOR-czasT %.2fsek.\n", czasSOR-czasT );

    std::cout << "END" << std::endl;
    return 0;
}