//
// Created by mateusz on 31.05.18.
//

#ifndef MO_11_PROJEKT_PROJEKT_H
#define MO_11_PROJEKT_PROJEKT_H
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <algorithm>
using namespace std;


class Projekt{

public:

Projekt(int ile_xx,bool zapis);
~Projekt(void){};

    void alokuj(double **tablica,int N, int M);
    void rozwiaz_analityczne();
    //double** rozwiaz_analityczne();
    void warunek(double **roz);
    void save_gnuplot(double **roz, string nazwa );
    //void saveA_gnuplot2( string nazwa, string roz );

    void rozwiaz_laasonen_thomasa();
    double** blad_bezwgl(double** blad, double **roz, double max_blad,string nazwa);
    //double MaksymalnyBladDlaDanegoT(double **U);

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


void rozwiaz_laasonen_SOR();


    void AlgorytmThomasa_rozwiazanie(double *nadDiag, double *Diag, double *podDiag, double *B, int n, double *X);

    void AlgorytmThomasa_macierz(double *nadDiag, double *Diag, double *podDiag, int n);

    void save_gnuplot2(double **roz,string nazwa, string rozsz);
    void save_gnuplot3(double **roz,string nazwa, string rozsz);
    void SOR(double **A, double *B, double *X0, int n);

    double max(double *x, int n, double **A, double *b);

    double Estymator(double *x0, double *x1, int n);
    void f_blad ();

    void zapis_do_pliku_w(char *nazwa, double *wektor);

    void saveA(double **roz, string nazwa);
};


#endif //MO_11_PROJEKT_PROJEKT_H
