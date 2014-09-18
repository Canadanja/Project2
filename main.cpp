// calculates the one-dimensional poisson equation
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <armadillo>

using namespace std;
using namespace arma;

double fillmatrix (int, double, double);
double searchmatrix (double, int, int, int);
void jacobi (double **, double **, int, int, int);


int main()
{
    unsigned int n, i;
    double x, x_min, x_max;
    double *u, u_min, u_max;
    double h_step, V;

    n=100;
    x_min = 0.;                                     // area of calculation    //
    x_max = 1.;                                     //                        //
    x = x_min;                                      // initial value for var. //
    u_min = 0.;                                     // dirichlet-boundaries   //
    u_max = 0.;                                     //                        //
    h_step = (x_max - x_min) / n;                   // stepsize               //
    V =                                             // potential              //


    fillmatrix(n, h_step, V);
    searchmatrix();
    jacobi();


    u = new double[n];                              // initialize array for   //
    u[0] = u_min;                                   // solution u             //
    u[n] = u_max;                                   //                        //

    // tridiag solves the equation Au = f
    tridiag_poisson (h_step, n, u);

    // output
    for (i = 0; i <= n + 1; i++)
    {
        cout << h_step*i << "    " << u[i] << endl;
    }

    delete [] u;
    return 0;
}

double fillmatrix(int n, double h_step, double V)
{
    mat A(n,n);
    A.zeros();
    A(0,0)=2/(h_step*h_step) + V;
    for(int i=1; i<n; i++){
        A(i,i) = 2/(h_step*h_step) + V;
        A(i,i-1) = -1/(h_step*h_step);
        A(i-1,i) = -1/(h_step*h_step);
    }
    return A;
}

double searchmatrix(double **A,int p, int q, int n)
{
    double max;
    for (int i = 0; i < n; ++i) {
        for ( int j = i+1; j < n; ++j) {
            double aij = fabs(A[i][j]);
            if ( aij > max)
                {
                max=aij; p=i;q=j;
                }
            }
        }
}

void jacobi(double **A, double **R, int k, int l, int n)
{
    ￼double s, c;
    if ( A[k][l] != 0.0 ) {
        double t, tau;
        tau = (A[l][l] −A[k][k])/(2∗A[k][l]);
        if ( tau >= 0 ) {
            t = 1.0/(tau + sqrt(1.0 + tau∗tau));
        }
        else {
            t = −1.0/(−tau +sqrt (1.0 + tau∗tau ) ) ;
        }
        c = 1/sqrt(1+t∗t);
        s = c∗t;
    }
    else {
        c=1.0;
        s=0.0;
    }

    double akk, all, aik, ail, rik, ril;
    akk = A[k][k];
    all =A[l][l];
    A[k][k] = c∗c∗a kk − 2.0∗c∗s∗A[k][l] + s∗s∗all ;
    A[l][l] = s∗s∗a kk + 2.0∗c∗s∗A[k][l] + c∗c∗all;
    A[k][l] = 0.0; // hard−coding non−diagonal elements by hand
    A[l][k] = 0.0; // same here
    for ( int i = 0; i < n; i++ ) {
        if ( i != k && i != l ) {
            aik = A[i][k];
            ail =A[i][l];
            A[i][k] = c∗aik − s∗ail;
            A[k][i] =A[i][k];
            A[i][l] = c∗ail + s∗aik;
            A[l][i] = A[i][l];
        }
        rik =R[i][k]; // new eigenvectors
        ril =R[i][l];
        R[i][k] = c∗rik − s∗ril;
        R[i][l] = c∗ril + s∗rik;
    }
    return ;
}

//double func (double h_step, double x)
//{
//    // second derivative rhs-function
//    double f;

//    f = h_step*h_step*100.*exp(-10.*x);
//    return f;
//}

// //TODO: Change the call of func and or func itself because of unnecessary FLOPs
//void tridiag (double h_step, int n, double *u){
//    /* Uses the Thomas-Algorithm, code from "Numerical Recipes,               *
//     * Third Edition", p. 56 f.                                               */
//    double a,b,c;
//    double btemp;
//    unsigned int i;
//    vec temp(n);

//    a = -1;
//    c = a ;
//    b = 2 ;

//    btemp = b;
//    u[1] = func(h_step, h_step*1.)/btemp;
//    for(i = 2 ; i <= n; i++){
//        temp[i] = c/btemp;
//        btemp = b - a*temp[i];
//        u[i] = (func(h_step, h_step*i) - a*u[i-1])/btemp;
//    }

//    for(i = n-1; i >= 1; i--){
//        u[i] -= temp[i+1]*u[i+1];
//    }
//}


//TODO: Change the call of func and or func itself because of unnecessary FLOPs
//void tridiag_poisson (double h_step, int n, double *u){
//    /* Just usable for poisson equation tridiagonal matrices. Reduces number  *
//     * of FLOPs to 6N.                                                        */
//    unsigned int i;
//    vec ftemp(n);

//    for (i = 2; i <= n; i++){
//        ftemp[i] = ftemp[i-1] + i*func(h_step, h_step*i);
//    }

//    for (i = n-1; i >= 1; i--){
//        u[i] = (ftemp[i] + u[i+1]*i)/(1. + i);
//    }
//}
