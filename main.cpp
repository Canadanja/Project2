// calculates the one-dimensional poisson equation
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <armadillo>
#include <iomanip>

using namespace std;
using namespace arma;

mat fillmatrix (unsigned int, double, vec);
double offnorm (mat);
void searchlargest (mat, unsigned int &, unsigned int&);
void jacobi (mat &, mat &, unsigned int, unsigned int);
void sorting (mat, mat, mat &, vec &);


int main()
{
    unsigned int n = 100.;                          // array size             //
    unsigned int p = 0.;                            //                        //
    unsigned int q = 0.;                            //                        //
    double rho_min = 0.;                            // min. calc. area        //
    double rho_max = 5.;                            // max. calc. area        //
    double rho = rho_min;                           //                        //
    double h_step = (rho_max-rho_min)/n;            // stepwidth              //
    double tolerance = 1e-10;                       //                        //
    double omega = 0.5;                             // strength of osc. pot.  //
    double offA;                                    //                        //
    vec V = zeros<vec>(n+1);                        // potential              // 
    mat A, R;                                       // Mainmat., Eigenvec.    // 
    mat T;                                          // sorted Eigenvectors    //
    vec lambda;                                     // sorted Eigenvalues     //
    mat D, eigvec;                                  // same for armadillo     //
    vec eigval;                                     //                        //  

    // define potential
    // TODO: write an init_potential function for arbitrary potentials
    for (unsigned int i = 0; i <= n; i++)
    {
        V(i) = rho*rho;
        //V(i) = rho*rho*omega*omega-1./rho;
        rho += h_step;
    }

   /*--------------------------------------------------------------------------/
   /                  Jacobi algorithm -- main part                            /          
   /--------------------------------------------------------------------------*/
    A = fillmatrix(n,h_step,V);                     // dense,symmetric matrix //
    R.eye(n-1,n-1);                                 // eigenvector-matrix     //

    offA = offnorm(A);
    while (offA > tolerance)
    {
        searchlargest(A,p,q);                       // search f. larg. elem.  // 
        jacobi(A,R,p,q);                            // jacobi rotation        // 
        offA = offnorm(A);                          // calculate new norm     // 
    }
    
    // bring the eigenvalues and its associated eigenvectors in right order   // 
    sorting(A,R,T,lambda);                          // T, lambda -> sorted    // 
    // TODO: Normalization is missing!
    /*------------------------------------------------------------------------*/

    // Armadillo routine for eigenvalues
    D = fillmatrix(n,h_step,V);                      
    eig_sym(eigval,eigvec,D);

    // output
    T.print();
//    cout << "Jacobi" << setw(10) << "Armadillo" << endl;
//    for (unsigned int i = 0; i < 9; i++)
//    {
//        cout << G(i) << setw(10) << eigval(i) << endl;
//    }

    return 0;
}

mat fillmatrix(unsigned int n, double h_step, vec V)
{
    double h_step2 = h_step*h_step;
    mat A = zeros<mat>(n-1,n-1);

    A(0,0)=2./h_step2 + V(1);
    for(unsigned int i = 1; i < n-1; i++){
        A(i,i) = 2./h_step2 + V(i+1);
        A(i,i-1) = -1./h_step2;
        A(i-1,i) = -1./h_step2;
    }
    return A;
}

double offnorm(mat A)
{
    unsigned int n = A.n_rows;                      // arraysize              //
    double offA;                                    // norm over nondiags     //

    // loop over (upper) nondiagonal elements
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = i+1; j < n; j++) {
                offA += 2.*A(i,j)*A(i,j);           // 2x for symmetric mat.  // 
        }
    }
    offA = sqrt(offA);

    return offA;
}

void searchlargest(mat A, unsigned int &p, unsigned int &q)
{
    unsigned int n = A.n_rows;                      // arraysize              //
    double max = 0.0;                               // largest value in mat.  //

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
    unsigned int n = A.n_rows;                      // arraysize              //
    double s, c;                                    // sin(theta),cos(theta)  // 

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
    all = A(l,l);
    A(k,k) = c*c*akk - 2.0*c*s*A(k,l) + s*s*all ;
    A(l,l) = s*s*akk + 2.0*c*s*A(k,l) + c*c*all;
    A(k,l) = 0.0; // hard-coding non-diagonal elements by hand
    A(l,k) = 0.0; // same here
    for (unsigned int i = 0; i < n; i++ ) {
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

void sorting(mat A, mat R, mat &T, vec &lambda)
{
    unsigned int n = A.n_rows;                      // assume symmetry        //

    // Sorting of eigenvalues
    vec K;
    K.zeros(n);
    for (unsigned int i = 0; i < n; i++)
    {
        K(i) = A(i,i);
    }
    lambda = sort(K);

    // Sorting of eigenvectors
    T.zeros(n,n);
    uvec ev_indices = sort_index(K);
    for (unsigned int i = 0; i < n; i++)
    {
        for (unsigned int j = 0; j < n; j++)
        {
            T(j,i) = R(j,ev_indices(i));
        }
    }
}
