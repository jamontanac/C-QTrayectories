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

// ------------------------------------------------------------

int main(int argc, char **argv)
{
  int Trajectories = 500;
  // int Trajectories=atoi(argv[1]);
  std::cout.precision(16);
  std::cout.setf(std::ios::scientific);
  std::mt19937 generator(1);
  std::normal_distribution<double> Noise(Mu, Sigma);//set a normal distribution to call it we do Noise(generator)
  int Size = 2;

  // -------------- Initial Conditions --------------
  // std::complex<double> x0 (1.0,0.0);
  // std::complex<double> y0 (0.0,0.0);
  // double t=0.0;

  // -------------------------------------------------

  // -----------------things to use----------------------
  arma::cx_mat rho1 (Size,Size);
  arma::cx_mat rho2 (Size,Size);
  arma::cx_mat ketstate (Size,1);
  arma::cx_mat dHamiltonian (Size,Size);
  arma::cx_mat drho (Size,Size);
  arma::vec X_avg(N);
  X_avg.fill(0.0);
  arma::vec Y_avg(N);
  Y_avg.fill(0.0);
  arma::vec Z_avg(N);
  Z_avg.fill(0.0);
  double dXi = 0.0;
  double Work = 0.0;
  double Energy = 0.0;
  double X_block=0.0;
  double Y_block=0.0;
  double Z_block=0.0;

  for(int trajectory = 0 ; trajectory< Trajectories;trajectory++ )
  {
    // -------------- Initial Conditions --------------
    std::complex<double> x0 (1.0,0.0);
    std::complex<double> y0 (0.0,0.0);
    double t=0.0;
    // std::cout<<"Start"<<"  "<<x0<<"  "<<y0<< "  "<< t<<std::endl;

// running a trajectory
    for(int step=0; step< N; step++)
      {
        rho1 = State(x0,y0);
        dXi = Noise(generator);
        ketstate = Euler(x0,y0,t,dt,dXi);
        ketstate = Normalise(ketstate);
        rho2 = State(ketstate(0,0),ketstate(1,0));
        dHamiltonian = Hamiltonian(t+dt)-Hamiltonian(t);
        drho = rho2 - rho1;
        Work = real(arma::trace(rho2*dHamiltonian));
        Energy = real(arma::trace(Hamiltonian(t+dt)*rho2)-arma::trace(Hamiltonian(t)*rho1));
        X_block = real(arma::trace(rho1*Sigma_x));
        X_avg[step] += X_block/Trajectories;
        Y_block = real(arma::trace(rho1*Sigma_y));
        Y_avg[step] += Y_block/Trajectories;
        Z_block = real(arma::trace(rho1*Sigma_z));
        Z_avg[step] += Z_block/Trajectories;
        // std::cout<<t/Tmax<<"  "<< X_block<<"  "<< Y_block<<"  "<< Z_block<<std::endl;
        t += dt;
        x0 = ketstate(0,0);
        y0 = ketstate(1,0);
      }
      // -------------- Initial Conditions --------------
      // x0= (1.0,0.0);
      // y0= (0.0,0.0);
      // t=0.0;
      // std::cout<<"end"<<"  "<<x0<<"  "<<y0<<"  "<<t<<std::endl;
    }
    double time=0.0;
    for(int i = 0;i<N;i++)
    {

      std::cout<<time/Tmax<<"  "<<X_avg[i]<<"  "<<Y_avg[i]<<"  "<<Z_avg[i]<<std::endl;
      time +=dt;
    }
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
  x1 += eulermethod(0,0);
  x2 += eulermethod(1,0);
  // return eulermethod;
    return Ket(x1,x2);
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
