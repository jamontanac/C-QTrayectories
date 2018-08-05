#include <iostream>
#include <cmath>
#include <complex>
#include "armadillo"
#include <vector>
// #include "omp.h"
const double Sigma=0.1;
const double Mu=0.0;
const double dt = 0.01;
const arma::cx_mat Sigma_y={{0.0,std::complex<double> (0.0,1.0)},{std::complex<double> (0.0,-1.0),0.0}};
const double Gamma = 5e-3;
std::complex<double> Correlated_noise (double Noise1, double Noise2);
int main(void)
{
  int N=10;
  std::mt19937 generator(1);
  std::normal_distribution<double> Noise1(Mu, Sigma);
  std::normal_distribution<double> Noise2(Mu, Sigma);
  for(int i=0; i<N; i++)
    {
      std::cout<<Correlated_noise(Noise1(generator),Noise2(generator))<<std::endl;
      
    }
  // std::cout<<inv(Eigenvectors)*D*Eigenvectors<<std::endl; //this is the way of making a diagonal matriz with the eigenvalues on its diagonal

  return 0;
}

std::complex<double> Correlated_noise (double Noise1, double Noise2)
{
  arma::vec Eigenvalues;
  arma::mat Eigenvectors;
  std::complex<double> U_corr (1.0,0.0);
  arma::vec target_mean  = {0.0, 0.0};
  arma::mat Matrix_Variance = {{1.0+real(U_corr),imag(U_corr)},{imag(U_corr),1.0-real(U_corr)}};
  Matrix_Variance=(dt/2)*Matrix_Variance;
  eig_sym(Eigenvalues, Eigenvectors, Matrix_Variance);
  arma::mat l = sqrt(diagmat(Eigenvalues));
  arma::mat Q = Eigenvectors * l;
  arma::vec original (2);
  original.zeros();
  arma::vec tweaked (2);
  original={Noise1,Noise2};
  // // std::cout << Q<<std::endl;
  tweaked= (Q*original + target_mean);
  std::complex<double> Cor_Noise (tweaked(0),tweaked(1));
  return Cor_Noise;
}
