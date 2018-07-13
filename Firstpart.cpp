#include <iostream>
#include <cmath>
#include <complex>
#include "armadillo"
#include "omp.h"
// --------------------- defining constants --------------------
const double Hbar = 1.0;
const double Epsilon=1.0;
const double Nu = 8.0;
const double Tau = 3.0;
const double Gamma = 5e-3;
const double G = 6.25;
const double dt = 0.01;
const double Tmax=Tau;
const int N=Tmax/dt;
const double Sigma=sqrt(dt);
const double Mu=0.0;
// ----------- Defining as constants the pauli matrices ----------------
const arma::cx_mat Sigma_x={{0.0,std::complex<double> (1.0,0.0)},{std::complex<double> (1.0,0.0),0.0}};
const arma::cx_mat Sigma_y={{0.0,std::complex<double> (0.0,1.0)},{std::complex<double> (0.0,-1.0),0.0}};
const arma::cx_mat Sigma_z={{std::complex<double> (-1.0,0.0), 0.0},{ 0.0,std::complex<double> (1.0,0.0)}};
const arma::cx_mat Sigma_plus={{0.0, 0.0},{ std::complex<double> (1.0,0.0),0.0}};
const arma::cx_mat Sigma_minous={{ 0.0,std::complex<double> (1.0,0.0)},{ 0.0,0.0}};
const arma::cx_mat J= sqrt(Gamma)*Sigma_z;
const arma::cx_mat Jdagga= J.t();
// ----------------------- Defining Functions ----------------------
double Lambda(double time);
// arma::cx_mat Hamiltonian(arma::cx_mat & H,double time);
arma::cx_mat Hamiltonian(double time);
arma::cx_mat Ket(std::complex<double> x1 , std::complex<double> x2);
arma::cx_mat Bra(std::complex<double> x1 , std::complex<double> x2);
arma::cx_mat Euler(std::complex<double> & x1,std::complex<double> &x2,double time,double dTime,double dW);
int main(int argc, char **argv)
{
  // int Trajectories=atoi(argv[1]);
  // arma::cx_mat J= sqrt(Gamma)*Sigma_y;
  // arma::cx_mat Jdagga=J.t();
  // std::cout<< J <<std::endl;
  // std::cout<< Jdagga <<std::endl;
  // std::cout<< Lambda(10)<<std::endl;
  // std::cout<<real(arma::trace(Sigma_x))<<std::endl;
  // std::cout<<real(arma::trace(Sigma_y))<<std::endl;
  // std::cout<<real(arma::trace(Sigma_z))<<std::endl;
  // std::cout<<Sigma_plus+Sigma_minous<<std::endl;
  // std::cout<<Sigma_minous<<std::endl;
  arma::mat A={{2,1},{1,1}};
  arma::mat B={{1,1},{3,1}};
  arma::mat C = A*B;
  std::cout<<C<<std::endl;
  // std::cout<<A<<std::endl;
  // A.fill(1.0);
  // std::cout<<A<<std::endl;

  // std::cout<<A*B<<std::endl;
  // std::cout<<real(Numero + Dos)<<std::endl;
  // std::cout<<imag(Numero + Dos)<<std::endl;
  
  // arma::cx_mat A(2,2);
  // std::cout<<A<<std::endl;
  // Hamiltonian(A,10);
  // std::cout<<A<<std::endl;
  // arma::cx_mat a=Ket(Numero,Dos);
  std::complex<double> Numero (1.0,2.0);
  std::complex<double> Dos (1.0,3.0);
  // arma:: cx_mat A=Euler(Numero,Dos,10,0.1,0.02);
  
  // arma::cx_mat A(2,2);
  // A.eye();
  // std::cout<<A<<std::endl;
  // std::cout<<Bra(Numero,Dos)<<std::endl;
  return 0;
}
// --------------- Implementing functions ----------------
double Lambda(double time)
{
  return G*(1.0/cosh(Nu*(1.0 - (time/Tau))));
}
arma::cx_mat Hamiltonian(double time)
{
  arma::cx_mat H = Epsilon*Sigma_z + Lambda(time)*Sigma_x;
  return H;
}
// arma::cx_mat Hamiltonian(arma::cx_mat & H,double time)
// {
//   H = Epsilon*Sigma_z + Lambda(time)*Sigma_x;
//   return H;
// }
arma::cx_mat Ket(std::complex<double> x1 , std::complex<double> x2)
{
  arma::cx_mat Matriz (2,1);
  Matriz(0,0) = x1;
  Matriz(1,0) = x2;
  return Matriz;
}
arma::cx_mat Bra(std::complex<double> x1 , std::complex<double> x2)
{
  arma::cx_mat Matriz (1,2);
  Matriz(0,0) = std::conj(x1);
  Matriz(0,1) = std::conj(x2);
  return Matriz;
}
arma::cx_mat Euler(std::complex<double> & x1,std::complex<double> & x2,double time,double dTime,double dW)
{
  arma::cx_mat Id(2,2);
  Id.eye();
  arma::cx_mat A1=(Id-Jdagga);
  return A1;
    
}
