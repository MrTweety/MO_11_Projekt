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

double *xxxT, *bmaxxxT;
double *xxxS,*bmaxxxS;


int main() {




    //*************************************************//
    //*******obliczanie max wartość bezwzglednaej******//
    //******dla t_max w funkcji kroku przestrz. h******//
    //*************************************************//


    std::cout << "Hello, World!" << std::endl;
    int count =6;
    int i=count;

    fstream bmaxT,bmaxTlog,bmaxS,bmaxSlog;
    bmaxT.open("W1_bmax_laasonen_thomasa.csv", ios::out);
    bmaxTlog.open("W1_bmax_laasonen_thomasa_log.csv", ios::out);
    bmaxS.open("W1_bmax_laasonen_SOR.csv", ios::out);
    bmaxSlog.open("W1_bmax_laasonen_SOR_log.csv", ios::out);

    for ( i =0 ; i < count ; i++) {

        cout<<"i==="<<i<<endl;
        int ile_iteracji = pow(2.0, i + 4);

        xxxT = new double[count];
        bmaxxxT = new double[count];
        xxxS = new double[count];
        bmaxxxS = new double[count];

        Projekt proj(ile_iteracji, 0);
        proj.rozwiaz_analityczne();
        proj.rozwiaz_laasonen_thomasa();
        proj.rozwiaz_laasonen_SOR();
        //proj.save_gnuplot3(proj.rozwiazanieT, "a", "a");
        //proj.f_blad();

        proj.maxb(proj.rozwiazanieT, xxxT, bmaxxxT, i);
        proj.maxb(proj.rozwiazanieSOR, xxxS, bmaxxxS, i);

        //proj.saveA(proj.blad_T, "bladT.csv");

        bmaxT<<xxxT[i]<<"\t"<<bmaxxxT[i]<<endl;
        bmaxTlog<<log10(fabs(xxxT[i]))<<"\t"<<log10(fabs(bmaxxxT[i]))<<endl;

        bmaxS<<xxxT[i]<<"\t"<<bmaxxxT[i]<<endl;
        bmaxSlog<<log10(fabs(xxxT[i]))<<"\t"<<log10(fabs(bmaxxxT[i]))<<endl;

    }

    bmaxT.close();
    bmaxTlog.close();
    bmaxS.close();
    bmaxSlog.close();




    //*************************************************//
    //*******obliczanie czasu działąnia algorytmu******//
    //*************************************************//

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