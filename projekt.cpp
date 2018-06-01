//
// Created by mateusz on 31.05.18.
//

#include "projekt.h"


Projekt::Projekt(int ile_xx,bool zapis)
{
    D=1.0;
    t_max=2.0;
    alfa=1.0;
    beta=0.0;
    // 1.
    x_min =-6.*sqrt(D*t_max)-1.0;
    x_max = -x_min;
    this->zapis=zapis;


    ile_x  = ile_xx;

    // 2. obliczenie kroku h i dt
    h  =  (x_max -x_min)/ (ile_x-1);

    dt = h * h;
    //dt = t_max/(ile_t -1);


    ile_t =(int)( t_max/dt +1);

    dt = (t_max/(ile_t-1));

    lambda = D * dt / (h*h);

    cout<<"ile_x = "<<ile_x<<"\tile_t = "<<ile_t<<"\nt_max = "<<t_max<<endl;
    cout<<"h = "<<h<<"\tkrok dt = "<<dt<<endl;
    cout<<"lambda = "<<lambda<<endl;


    rozwiazanieT  = new double*[ile_t];
    for (int i=0;i<ile_t;i++)
        rozwiazanieT[i] = new double[ile_x];


    rozwiazanieA  = new double*[ile_t];
    for (int i=0;i<ile_t;i++)
        rozwiazanieA[i] = new double[ile_x];

    rozwiazanieSOR = new double*[ile_t];
    for (int i=0;i<ile_t;i++)
        rozwiazanieSOR[i] = new double[ile_x];


}


void Projekt::warunek(double **roz)
{
    double x ;
    for( int i = 0; i < ile_x; i++ ){ //warunek poczatkowy

        x=h*(double)i+x_min;


        if (x>=0)
            roz[0][i] = 0.0;

        else if   (x<0)
            roz[0][i] = 1.0;




    }

    for( int i = 0; i < ile_t; i++ ) {// warunki brzegowe , dla kzdej chwili T - konkretna wartosc 1. i ostatniej zmiennej X
        roz[i][0] = alfa;
        roz[i][ile_x-1] = beta;
    }

}




//double** Projekt::rozwiaz_analityczne()
void Projekt::rozwiaz_analityczne()
{
    double x,t;

    warunek(rozwiazanieA);

    //cout<<rozwiazanieA;


    x=x_min+h;
    t=dt;

    for (int i=1;i<ile_t;i++)
    {
        for (int j=1 ; j<ile_x-1;j++)
        {

            rozwiazanieA[i][j] = 0.5 * erfc( x / (2.0 * sqrt(D*t)));

            x+=h;
        }
        x=x_min+h;
        t+=dt;


    }
    if(zapis) {
        save_gnuplot(rozwiazanieA, "rozwA.txt");
        save_gnuplot2(rozwiazanieA, "rozwA", ".txt");
    }
}





void Projekt::rozwiaz_laasonen_thomasa(){
    cout<<"rozwiaz_laasonen_thomasa():"<<endl;
    double *nadDiag = new double[ile_x-1];
    double *Diag = new double[ile_x];
    double *podDiag = new double[ile_x-1];
    double *b = new double[ile_x];
    double *wyn = new double[ile_x];


    warunek(rozwiazanieT);

    nadDiag[0]  =lambda;
    Diag[0]=1.0;
    podDiag[0]=lambda;



    for(int i=1 ; i<ile_x-1; i++){
        nadDiag[i]   = lambda;
        Diag[i] =-( 1 + (2*lambda) );
        podDiag[i] = lambda;

    }

    Diag[ile_x-1]=1.0;


    //pierwszy etap alg Thomassa - zerowanie podprzekatnej:

    AlgorytmThomasa_macierz( nadDiag,  Diag,  podDiag, ile_x);


    // przebieg czasowy:
    for( int k = 1; k < ile_t; k++ ) {

        b[0] = alfa; //wynika z lewego warunku brzegowego

        for(int i=1 ; i<ile_x; i++)
            b[i] = - rozwiazanieT[k-1][i];

        b[ile_x-1] = beta ;

        AlgorytmThomasa_rozwiazanie(nadDiag,  Diag, podDiag, b, ile_x, wyn);


        //kopiowanie wynikow z chwili "k'atej" do macierzy wynikow:
        for( int i = 1; i < ile_x-1; i++ )
            rozwiazanieT[k][i] = wyn[i];


    }
    if(zapis) {
        save_gnuplot(rozwiazanieT, "rozwT.txt");
        save_gnuplot2(rozwiazanieT, "rozwT", ".txt");
    }
}


void Projekt::AlgorytmThomasa_macierz(double *nadDiag, double *Diag, double *podDiag, int n) {

    double l;
    int i;

    for (i = 1; i < n; i++)
    {
        //l = podDiag[i - 1] / Diag[i - 1];
        Diag[i] = Diag[i] - (podDiag[i - 1] / Diag[i - 1]) * nadDiag[i - 1];//wyliczamy eta[i]

    }
}

void Projekt::AlgorytmThomasa_rozwiazanie(double *nadDiag, double *Diag, double *podDiag, double *B, int n, double *X) {

    double l;
    int i;

    for (i = 1; i < n; i++)
    {
        //l = podDiag[i - 1] / Diag[i - 1];
        B[i] = B[i] - (podDiag[i - 1] / Diag[i - 1])* B[i - 1];//wyliczamy r[i]
    }

    X[n-1] = B[n-1] / Diag[n-1];//trójkatna górna
    for (i = n - 2; i >= 0; i--)
        X[i] = (B[i] - nadDiag[i] * X[i + 1]) / Diag[i]; //wyliczamy X[i]

}

void Projekt::rozwiaz_laasonen_SOR(){
    double x;
    double *b = new double[ile_x];
   // double *wyn = new double[ile_x];
    double *X0 = new double[ile_x];//przyblizenie pcozatkowe [0,0,...0,0]



    x = x_min;
    for( int i = 0; i < ile_x; i++ ){ //warunek poczatkowy

        if(x<0)
            rozwiazanieSOR[0][i] = 1.0;
        else
            rozwiazanieSOR[0][i] = 0.0;
        x=x+h;
    }


    for( int i = 0; i < ile_t; i++ ) {// warunki brzegowe
        rozwiazanieSOR[i][0] = alfa;
        rozwiazanieSOR[i][ile_x-1] = beta;
    }

    double **A= new double*[ile_x];
    for (int i=0;i<ile_x;i++)
        A[i] = new double[ile_x];//macierz trojdiagonalna

//    for (int i=0;i<ile_x;i++) {
//        X0[i] = 0.0;
//        for (int j = 0; i < ile_x; j++)
//            A[i][j] = 0.0;
//    }


    A[1][0]=lambda;
    A[0][1]=lambda;
    A[0][0]=1.0;
    for( int i = 1; i < ile_x-1; i++ )//wypelnianie macierzy
    {
        A[i+1][i] = lambda;//down
        A[i][i] = -( 1 + (2*lambda) );//diag
        A[i][i+1] = lambda;//up
    }
    A[ile_x-1][ile_x-1] = 1.0;

    //przebieg czasowy
    for( int k = 1; k < ile_t; k++ )
    {


        b[0] =  alfa ;
        for( int i = 1; i < ile_x-1; i++ )//wypelnianie macierzy i
            b[i] = -rozwiazanieSOR[k-1][i];
        b[ile_x-1] = beta ;

        SOR(A,b,X0,ile_x);

        for( int i = 1; i < ile_x-1; i++ )//kopiowanie wynikow
            rozwiazanieSOR[k][i] = X0[i];
    }

    if(zapis) {
        save_gnuplot(rozwiazanieSOR, "rozwSOR.txt");
        save_gnuplot2(rozwiazanieSOR, "rozwSOR", ".txt");
    }
}


void Projekt::SOR(double **A, double *B, double *X0, int n) {
    //cout << "\n\nMetoda SOR\n\n";
    //cout << "i \t x[0] \t x[1] \t x[2] \t x[3] \t\t estymator \t\t residuum\n" << endl;

    double *Xp = new double[n];
    double estymator = 1, residuum = 1;
    double suma = 0;
    double Xtmp = 0;

    double omega = 0.5;
    int obieg = 0;

    double TOLX = 1e-9;
    double TOLF = 1e-9;
    int NMAX = 50;


    while ((residuum > TOLF && estymator > TOLX) && NMAX > obieg) {
        obieg++;

        for (int i = 0; i < n; i++) {
            suma = 0;

            for (int j = 0; j < n; j++)
                if (j != i)
                    suma += A[i][j] * X0[j];

            Xtmp= (B[i] - suma) / A[i][i];

            Xp[i] = X0[i];
            X0[i] = (1-omega)*X0[i] +omega*Xtmp;
        }

        estymator = Estymator(Xp, X0, n);
        residuum = max(X0, n, A, B);
    }
}



double Projekt::max(double *x, int n, double **A,double *b) {


    double *Res = new double[n];//r=||Ax-b||


    for (int i = 0; i < n; ++i) {
        Res[i] = 0;
        for (int j = 0; j < n; ++j)
            Res[i] += A[i][j] * x[j];
        Res[i] -= b[i];
        Res[i] = fabs(Res[i]);
    }



    double max = Res[0];
    for (int i = 0; i < n; ++i)
        if (fabs(Res[i]) > max)
            max = fabs(Res[i]);

    return max;
}
//********************************************************************
double Projekt::Estymator(double *x0, double *x1, int n) {

    for (int i = 0; i < n; ++i)
        x0[i]= fabs(x0[i] - x1[i]);;

    double max = x0[0];
    for (int i = 0; i < n; ++i)
        if (x0[i] > max)
            max = fabs(x0[i]);

    return max;
}










































void Projekt::save_gnuplot( double **roz,string nazwa )
{


    fstream fl;
    fl.open(nazwa.c_str(), ios::out);
    if (!fl.is_open())
    {
        cout<<"Blad otwierania pliku!"<<endl;
        system("pause");
        return;
    }


    double t,x;
    //w formie macierzy:



    x=x_min;
    t=0.0;
    for (int i=0;i<ile_t;i++,t=dt*i)
    {


        for(int j=0 ;  j<ile_x ; j++,x=x_min + h*j ){

            fl<<t<< "\t" <<  x  <<"\t"<<roz[i][j]<<endl;
        }
        fl<<endl;
        x=x_min;

    }



    fl.close();
    cout<<"Zapisano do pliku :  \""<<nazwa.c_str()<<"\""<<endl;

}



void Projekt::save_gnuplot2(double **roz, string nazwa, string rozsz )
{

    double t,x;
    //w formie macierzy:



    x=x_min;
    t=0.0;
    for (int i=0;i<ile_t;i++,t=dt*i)
    {
        if(i==0|| i==50 || i==100|| i==200|| i==300 || i==400 || i==511) {

            string nazwa2 = nazwa + to_string(i) + rozsz;
            fstream fl;
            fl.open(nazwa2.c_str(), ios::out);
            if (!fl.is_open()) {
                cout << "Blad otwierania pliku!" << endl;
                system("pause");
                return;
            }


            for (int j = 0; j < ile_x; j++, x = x_min + h * j) {

                fl << t << "\t" << x << "\t" << roz[i][j] << endl;
            }
            fl << endl;
            x = x_min;

            fl.close();
            cout << "Zapisano do pliku :  \"" << nazwa2.c_str() << "\"" << endl;
        }
    }

}


