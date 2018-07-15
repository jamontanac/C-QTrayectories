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
arma::cx_mat Hamiltonian(double time);
arma::cx_mat Ket(std::complex<double> x1 , std::complex<double> x2);
arma::cx_mat Bra(std::complex<double> x1 , std::complex<double> x2);
arma::cx_mat Euler(std::complex<double> & x1,std::complex<double> &x2,double time,double dTime,double dW);
arma::cx_mat State(std::complex<double> x1,std::complex<double> x2);
arma::cx_mat Normalise(arma::cx_mat Not_normalised);
int main(int argc, char **argv)
{
  std::cout.precision(16);
  std::cout.setf(std::ios::scientific);
  std::mt19937 generator(1);
  std::normal_distribution<double> Noise(Mu, Sigma);//set a normal distribution to call it we do Noise(generator)
  int Size = 2;
  
  // -------------- Initial Conditions --------------
  std::complex<double> x0 (1.0,0.0);
  std::complex<double> y0 (0.0,0.0);
  double t=0.0;

  // -------------------------------------------------

  // -----------------thing to use----------------------
  arma::cx_mat rho1 (Size,Size);
  arma::cx_mat rho2 (Size,Size);
  arma::cx_mat ketstate (Size,1);
  double dXi = 0.0;
  
  // int Trajectories=atoi(argv[1]);

  // ---------------------------------------------------------
  // std::mt19937 generator(std::random_device{}()); //std::mersenne_twister_engine set with seed the time of the computer with _64 we set to 64 bits
  // std::uniform_int_distribution<int> distribution(-100,100);//make a uniform distribution
  // ---------------------------------------------------------------

   rho1= State(x0,y0);
   dXi=Noise(generator);
   ketstate=Euler(x0,y0,t,dt,dXi);
   std::cout<<ketstate<<std::endl;
   std::cout<<ketstate.t()*ketstate<<std::endl;
   ketstate=Normalise(ketstate);
   std::cout<<ketstate<<std::endl;
   std::cout<<ketstate.t()*ketstate<<std::endl;
   // std::cout<<rho1<<std::endl;


  
  // std::complex<double> Numero (1.0,2.0);
  // std::complex<double> Dos (1.0,3.0);
  
  // arma:: cx_mat A=Euler(Numero,Dos,10,0.1,0.02);
  // arma:: cx_mat A=State(Numero,Dos);
  // std::cout<<arma::as_scalar(Bra(Numero,Dos)*Ket(Numero,Dos))<<std::endl;
  // std::cout<<A(0,0)<<"  "<< A(1,0)<<std::endl;
  // std::cout<<Normalise(A)<<std::endl;
  // std::cout<<A<<std::endl;
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
  std::complex<double> Iimag(0.0,1.0);
  arma::cx_mat A = (Id - Jdagga*J*0.5 - (Iimag/Hbar)*Hamiltonian(time))*Ket(x1,x2);
  arma::cx_mat B = J*Ket(x1,x2);
  arma::cx_mat eulermethod = A*dTime+B*dW;
  return eulermethod;
}
arma::cx_mat State(std::complex<double> x1,std::complex<double> x2)
{
  arma::cx_mat rho = Ket(x1,x2)*Bra(x1,x2);
  return rho;  
}
arma::cx_mat Normalise(arma::cx_mat Not_normalised)
{
  std::complex<double> Norm=std::sqrt(arma::as_scalar(Not_normalised.t()*Not_normalised));
  arma::cx_mat normalisedstate=Not_normalised/Norm;
  return normalisedstate;
}
