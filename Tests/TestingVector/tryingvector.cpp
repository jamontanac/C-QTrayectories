#include <iostream>
#include <cmath>
#include <complex>
#include "armadillo"
#include "omp.h"
#include <vector>
int main(void)
{
  int N=400;
  arma::cx_vec a (N);
  a.fill(0.0);
  // std::cout<< a <<std::endl;
  for(int i=0;i<N;i++)
  {
    a[i]=std::complex<double> (i,i+1);
  }
  // std::cout<< a <<std::endl;
  std::vector < std::complex <double> >A (N);
  for(int i=0;i<N;i++)
  {
    A[i]=std::complex<double> (i,i+1);
    // std::cout<<A[i]<<"  ";
  }
  // arma::cx_mat C = arma::randu <arma::cx_mat> (3,3);
  arma::cx_mat D = {{1,std::complex<double> (2,3)},{std::complex<double> (2,-3),2}};
  // std::cout<<C<<std::endl;
  // std::cout<<D<<std::endl;
  arma::vec Eigenvalues;
  arma::cx_mat Eigenvectors;
  eig_sym(Eigenvalues,Eigenvectors, D);
  // std::cout<< D<<std::endl;
  std::cout<< Eigenvalues<<std::endl;
  std::cout<< Eigenvectors<<std::endl;
  // Eigenvectors.print("Eigen vectors");
  double Eigenvalue1=0.0;
  double Eigenvalue2=0.0;
  arma::cx_vec Eigenvector1 (2);
  arma::cx_vec Eigenvector2 (2);
  for(int i =0;i<2;i++)
  {
    Eigenvalue1=Eigenvalues[0];
    Eigenvalue2=Eigenvalues[1];
    Eigenvector1[i]=Eigenvectors(i,0);
    Eigenvector2[i]=Eigenvectors(i,1);
  }
  // Eigenvector1.print('primero');
  std::cout<<Eigenvector1<<std::endl;
  std::cout<<Eigenvector2<<std::endl;
  std::cout<<(D*Eigenvector1-Eigenvalue1*Eigenvector1)<<std::endl;
  std::cout<<(D*Eigenvector2-Eigenvalue2*Eigenvector2)<<std::endl;
  return 0;

}
