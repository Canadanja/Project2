// calculates the one-dimensional poisson equation
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <armadillo>

using namespace std;
using namespace arma;

mat fillmatrix (unsigned int, double, vec);
double offnorm (mat);
void searchlargest (mat, unsigned int &, unsigned int&);
void jacobi (mat &, mat &, unsigned int, unsigned int);


int main()
{
    unsigned int n, p, q;
    double rho_min, rho_max;
    double rho;
    double h_step, offA;
    double tolerance;
    mat A, R;                                       // Mainmat., Eigenvec.    // 
    vec V;
    

    n=10;                                          // dimension              //
    tolerance = 1e-10;                              // tolerance              //
    rho_min = 0.;                                   // area of calculation    //
    rho_max = 5.0;                                   //                        //
    rho = rho_min;
    h_step = (rho_max-rho_min)/n;
    p = 0;
    q = 0;     

    // define potential
    V.zeros(n+1);
    for (unsigned int i = 0; i <= n; i++)
    {
        V(i) = rho*rho;
        rho += h_step;
    }

   // Jacobi main
    A = fillmatrix(n, h_step, V);
    R.eye(n-1,n-1);
    offA = offnorm(A);
    while (offA > tolerance)
    {
        searchlargest(A,p,q);                       // search f. larg. elem.  // 
        jacobi(A,R,p,q);                            // jacobi rotation        // 
        offA = offnorm(A);                          // calculate new norm     // 
    }

    // output
//    for (unsigned int i = 0; i < n-1; i++)
//    {
//        cout << A(i,i) << endl;
//    }
    A.print();

    return 0;
}

mat fillmatrix(unsigned int n, double h_step, vec V)
{
    mat A(n-1,n-1);
    A.zeros(); 
    A(0,0)=2./(h_step*h_step) + V(1);
    for(unsigned int i = 1; i < n-1; i++){
        A(i,i) = 2./(h_step*h_step) + V(i+1);
        A(i,i-1) = -1./(h_step*h_step);
        A(i-1,i) = -1./(h_step*h_step);
    }
    return A;
}

double offnorm(mat A)
{
    unsigned int n;                                 // arraysize              //
    double offA;                                    // norm over nondiags     //

    n = A.n_rows;                                   // assume symmetry        //    

    // loop over (upper) nondiagonal elements
    //TODO: Get certain about indices. (Compare to searchlargest routine)
    for (unsigned int i = 1; i < n; i++) {
        for (unsigned int j = i+1; j < n; j++) {
                offA += 2.*A(i,j)*A(i,j);           // 2x for symmetric mat.  // 
        }
    }
    offA = sqrt(offA);

    return offA;
}

void searchlargest(mat A, unsigned int &p, unsigned int &q)
{
    unsigned int n;                                 // arraysize              //
    double max;                                     // largest value in mat.  //
                                                                                
    n  = A.n_rows;                                  // assume symmetry        //
    max = 0.;

    // loop over (upper) nondiagonal elements 
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = i+1; j < n; j++) {
            double aij = fabs(A(i,j));
            if ( aij > max)
                {
                max=aij; p=i;q=j;
                }
            }
        }
}

void jacobi(mat &A, mat &R, unsigned int k, unsigned int l)
{
    unsigned int n;                                 // arraysize              //
    double s, c;                                    // sin(theta),cos(theta)  // 

    n  = A.n_rows;                                  // assume symmetry        //

    if ( A(k,l) != 0.0 ) {
        double t, tau;                              // for defining s,c       //
        tau = (A(l,l) - A(k,k))/(2.*A(k,l));

        if ( tau > 0 ) {
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        } else {
            t = -1.0/(-tau + sqrt(1.0 + tau*tau));
        }

        c = 1. / sqrt(1.+t*t);
        s = c*t;
    }
    else {
        c=1.0;
        s=0.0;
    }

    double akk, all, aik, ail, rik, ril;
    akk = A(k,k);
    all =A(l,l);
    A(k,k) = c*c*akk - 2.0*c*s*A(k,l) + s*s*all ;
    A(l,l) = s*s*akk + 2.0*c*s*A(k,l) + c*c*all;
    A(k,l) = 0.0; // hard-coding non-diagonal elements by hand
    A(l,k) = 0.0; // same here
    for (unsigned int i = 1; i < n; i++ ) {
        if (i != k && i != l ) {
            aik = A(i,k);
            ail = A(i,l);
            A(i,k) = c*aik - s*ail;
            A(k,i) = A(i,k);
            A(i,l) = c*ail + s*aik;
            A(l,i) = A(i,l);
        }
        rik =R(i,k); // new eigenvectors
        ril =R(i,l);
        R(i,k) = c*rik - s*ril;
        R(i,l) = c*ril + s*rik;
    }
    return ;
}

