#include <iostream>
#include <cmath>
#include <complex>
#include "armadillo"
#include <vector>
// #include "omp.h"
const double Sigma=0.1;
const double Mu=0.0;
const arma::cx_mat Sigma_y={{0.0,std::complex<double> (0.0,1.0)},{std::complex<double> (0.0,-1.0),0.0}};
const double Gamma = 5e-3;
std::complex<double> Correlated_noise (double Noise1, double Noise2);
int main(void)
{
  int N=10;
  std::mt19937 generator(1);
  std::normal_distribution<double> Noise1(Mu, Sigma);
  std::normal_distribution<double> Noise2(Mu, Sigma);
  // // std::vector<double> X(N);
  // // std::vector<double> Y(N);
  // // arma::vec X_arma(N);
  // // arma::vec Y_arma(N);
  // arma::vec Eigenvalues;
  // arma::mat Eigenvectors;
  // arma::mat D ={{1.0,0.0},{0.0,0.1}};
  // eig_sym(Eigenvalues, Eigenvectors, D);
  // arma::mat l = sqrt(diagmat(Eigenvalues));
  // arma::mat Q = Eigenvectors * l;
  // arma::vec original (2);
  // arma::vec target_mean  = {1.0, 5.0};
  // arma::vec tweaked (2);
  for(int i=0; i<N; i++)
    {
      std::cout<<Correlated_noise(Noise1(generator),Noise2(generator))<<std::endl;
      // original={Noise1(generator),Noise2(generator)};
      // tweaked= (Q*original + target_mean);
      // std::cout<<tweaked(0)<<"\t"<<tweaked(1)<<std::endl;
      
    }
  // std::cout<<inv(Eigenvectors)*D*Eigenvectors<<std::endl; //this is the way of making a diagonal matriz with the eigenvalues on its diagonal

  // std::cout<< l <<std::endl;
  // std::cout<< Q <<std::endl;
  
  
  










  // for(int i =0; i<2;i++)
  //   {
  //     for(int j=0;j<2;j++)
  // 	{
  // 	  std::cout<<D(i,j)<<"\t"<<D[i,j]<<std::endl;
  // 	}
  //   }
  // std::cout<<Eigenvalues<<"\t"<<Eigenvectors<<std::endl;
  
  // double Eigenvalue1=0.0;
  // double Eigenvalue2=0.0;
  // arma::vec Eigenvector1 (2);
  // arma::cx_vec Eigenvector2 (2);
  // for(int i = 0 ; i< N; i++)
  //   {
  //     X[i] = Noise1(generator);
  //     Y[i] = Noise2(generator);
  //     X_arma[i] = Noise1(generator);
  //     Y_arma[i] = Noise2(generator);
  //     // std::cout<<Noise1(generator)<<"\t"<<Noise2(generator)<<std::endl;
  //   }
  // std::cout<<X_arma<<" otro"<<Y_arma<<std::endl;
  
  return 0;
}

std::complex<double> Correlated_noise (double Noise1, double Noise2)
{
  // std::mt19937 generator(1);
  // std::normal_distribution<double> Noise1(Mean, Variance);
  // std::normal_distribution<double> Noise2(Mean, Variance);
  arma::vec Eigenvalues;
  arma::mat Eigenvectors;
  arma::mat Matrix_Variance ={{1.0,0.7},{0.7,1.0}};
  arma::vec target_mean  = {0.0, 0.0};
  eig_sym(Eigenvalues, Eigenvectors, Matrix_Variance);
  arma::mat l = sqrt(diagmat(Eigenvalues));
  arma::mat Q = Eigenvectors * l;
  arma::vec original (2);
  original.zeros();
  arma::vec tweaked (2);
  original={Noise1,Noise2};
  // std::cout << Q<<std::endl;
  tweaked= (Q*original + target_mean);
  std::complex<double> Cor_Noise (tweaked(0),tweaked(1));
  return Cor_Noise;
  // std::cout<<tweaked(0)<<"\t"<<tweaked(1)<<std::endl;
}
