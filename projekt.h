//
// Created by mateusz on 31.05.18.
//

#ifndef MO_11_PROJEKT_PROJEKT_H
#define MO_11_PROJEKT_PROJEKT_H
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <algorithm> //replace
using namespace std;

class Projekt{

public:

    double **rozwiazanieT;
    double **rozwiazanieSOR;
    double **rozwiazanieA;

    double **blad_T;
    double **blad_SOR;
    double  blad_max_T;
    double  blad_max_SOR;

    int ile_t;
    int ile_x;

    double alfa;
    double beta;
    double lambda;
    double h;
    double dt;
    double D;
    double x_min;
    double x_max;
    double t_max;
    bool zapis;



    Projekt(int ile_xx,bool zapis);
    ~Projekt(void){};


    void rozwiaz_analityczne();
    void warunek(double **roz);
    void rozwiaz_laasonen_thomasa();
    void AlgorytmThomasa_macierz(double *nadDlag, double *Diag, double *podDiag, int n);
    void AlgorytmThomasa_rozwiazanie(double *nadDiag, double *Diag, double *podDiag, double *B, int n, double *X);
    void rozwiaz_laasonen_SOR();
    void SOR(double **A, double *B, double *X0, int n);
    double max(double *x, int n, double **A, double *b);
    double Estymator(double *x0, double *x1, int n);
    double** blad_bezwgl(double** blad, double **roz, double max_blad,string nazwa);
    void save_gnuplot_w3(double **roz,string nazwa, string rozsz);
    void maxb(double **roz,  double *tablica, double *tablica2, int i);
    void save_macierz(double **roz, string nazwa);

    void save_gnuplot_w2(double **roz, string nazwa, string rozsz);

    void save_gnuplot(double **roz, string nazwa);
};

#endif //MO_11_PROJEKT_PROJEKT_H
