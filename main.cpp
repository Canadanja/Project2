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
    double h_step, V;

    n=100;                                          // dimension              //
    tolerance = 1e-10;                              // tolerance              //
    rho_min = 0.;                                   // area of calculation    //
    rho_max = 1.;                                   //                        //
    rho = rho_min;                                  // initial value for var. //
    h_step = (x_max - x_min) / n;                   // stepsize               //
    V = rho**2.                                     // potential              //


    // Jacobi algorithm
    A = fillmatrix(n, h_step, V);
    offA = offnorm(A)

    while (offA > tolerance)
    {
        searchlargest(A,p,q);                            // search f. larg. elem. // 
        jacobi(A);                                  // jacobi rotation        // 
        offA = offnorm(A)                           // calculate new norm     // 
    }

    // output
    for (i = 0; i <= n + 1; i++)
    {
        cout << h_step*i << "    " << u[i] << endl;
    }

    delete [] u;
    return 0;
}

mat fillmatrix(int n, double h_step, double V)
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

mat offnorm(mat A)
{
    unsigned int n;                                 // arraysize              //
    double offA;                                    // norm over nondiags     //

    n = A.n_rows;                                   // assume symmetry        //    

    // loop over (upper) nondiagonal elements
    //TODO: Get certain about indices. (Compare to searchlargest routine)
    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
                offA += 2.*A(i,j)*A(i,j);           // 2x for symmetric mat.  // 
        }
    }
    offA = sqrt(offA)

    return offA
}

double searchlargest(mat A,int p, int q)
{
    unsigned int n;                                 // arraysize              //
    double max;                                     // largest value in mat.  //
                                                                                
    n  = A.n_rows;                                  // assume symmetry        //

    // loop over (upper) nondiagonal elements 
    for (int i = 0; i < n; ++i) {
        for (int j = i+1; j < n; ++j) {
            double aij = fabs(A(i,j));
            if ( aij > max)
                {
                max=aij; p=i;q=j;
                }
            }
        }
}

void jacobi(mat A, mat R, int k, int l)
{
    unsigned int n;                                 // arraysize              //
    double s, c;                                    // sin(theta),cos(theta)  // 

    n  = A.n_rows;                                  // assume symmetry        //

    if ( A(k,l) != 0.0 ) {
        double t, tau;                              // for defining s,c       //
        tau = (A(l,l) − A(k,k))/(2.∗A(k,l));
        if ( tau >= 0.0 ) {
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
    akk = A(k,k);
    all =A(l,l);
    A(k,k) = c∗c∗a kk − 2.0∗c∗s∗A(k,l) + s∗s∗all ;
    A(l,l) = s∗s∗a kk + 2.0∗c∗s∗A(k,l) + c∗c∗all;
    A(k,l) = 0.0; // hard−coding non−diagonal elements by hand
    A(l,k) = 0.0; // same here
    for ( int i = 0; i < n; i++ ) {
        if ( i != k && i != l ) {
            aik = A(i,k);
            ail =A(i,l);
            A(i,k) = c∗aik − s∗ail;
            A(k,i) =A(i,k);
            A(i,l) = c∗ail + s∗aik;
            A(l,i) = A(i,l);
        }
        rik =R(i,k); // new eigenvectors
        ril =R(i,l);
        R(i,k) = c∗rik − s∗ril;
        R(i,l) = c∗ril + s∗rik;
    }
    return ;
}
