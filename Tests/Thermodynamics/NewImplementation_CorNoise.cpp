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
// const double Sigma=sqrt(dt);
const double Sigma = 1.0;
const double Mu=0.0;

// -------------------- Set the size of our problem ----------------
const int Size = 2;
// -------------- Initial Conditions as constants --------------
const std::complex<double> X_initial (1.0,0.0);
const std::complex<double> Y_initial (0.0,0.0);

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
std::complex<double> Correlated_noise (double Noise1, double Noise2);
arma::cx_mat Hamiltonian(double time);
arma::cx_mat Ket(std::complex<double> x1 , std::complex<double> x2);
arma::cx_mat Bra(std::complex<double> x1 , std::complex<double> x2);
arma::cx_mat Euler(arma::cx_mat & VectorState,double time,double dTime, std::complex<double> dW);
arma::cx_mat Unitary(arma::cx_mat & VectorState,double time,double dTime, std::complex<double> dW);
arma::cx_mat Non_Unitary(arma::cx_mat & VectorState,double time,double dTime, std::complex<double> dW);
arma::cx_mat State(arma::cx_mat VectorState);
arma::cx_mat Normalise(arma::cx_mat Not_normalised);
arma::cx_mat Initial_Conditions(std::complex<double> x1,std::complex<double> x2);

// ------------------------------------------------------------

int main(int argc, char **argv)
{
  int Trajectories = 1;
  // int Trajectories=atoi(argv[1]);
  std::cout.precision(16);
  std::cout.setf(std::ios::scientific);
  std::mt19937 generator(1);
  std::normal_distribution<double> Noise1(Mu, Sigma);
  std::normal_distribution<double> Noise2(Mu, Sigma);
  // std::normal_distribution<double> Noise(Mu, Sigma);//set a normal distribution to call it we do Noise(generator)
  
  
  
  // -----------------things to use----------------------
  arma::cx_mat New_ket(Size,1);
  arma::cx_mat rho0 (Size,Size);
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
  std::complex<double> dXi (0.0,0.0);
  double Work = 0.0;
  double Energy = 0.0;
  double Heat = 0.0;
  double X_block=0.0;
  double Y_block=0.0;
  double Z_block=0.0;
  // -------------------- Start the method to run over trajetories -----------
  for(int trajectory = 0 ; trajectory< Trajectories;trajectory++ )
    {
      
      New_ket = Initial_Conditions(X_initial,Y_initial);
      double t = 0.0;
      
      // running a trajectory
      for(int step=0; step< N; step++)
	{
	  rho0 = State(New_ket);
	  dXi=Correlated_noise(Noise1(generator),Noise2(generator));
	  // dXi = Noise(generator);
	  // --------------- Old Method --------------------
	  // ketstate = Euler(New_ket,t,dt,dXi);
	  // ----------- Evolve unitarily ----------------
	  ketstate = Unitary(New_ket,t,dt,dXi);
	  ketstate = Normalise(ketstate);
	  rho1 = State(ketstate);
	  // ------------------Evolve non unitarily ---------------
	  ketstate = Non_Unitary(ketstate,t,dt,dXi);
	  ketstate = Normalise(ketstate);
	  rho2= State(ketstate);

	  dHamiltonian = Hamiltonian(t+dt)-Hamiltonian(t);
	  drho = rho2 - rho0;
	  
	  Work = real(arma::trace(rho0*dHamiltonian));
	  Heat = real(arma::trace((rho2-rho1)*Hamiltonian(t)));
	  Energy = real(arma::trace( (dHamiltonian*rho2))) +real(arma::trace(drho*Hamiltonian(t)));
	  
	  X_block = real(arma::trace(rho0*Sigma_x));
	  X_avg[step] += X_block/Trajectories;
	  Y_block = real(arma::trace(rho0*Sigma_y));
	  Y_avg[step] += Y_block/Trajectories;
	  Z_block = real(arma::trace(rho0*Sigma_z));
	  Z_avg[step] += Z_block/Trajectories;
	
	  // std::cout<<t/Tmax<<"  "<< X_block<<"  "<< Y_block<<"  "<< Z_block<<std::endl;
	  std::cout<<t/Tmax<<"  "<< Energy<<"  "<< Heat <<"  "<< Work<<std::endl;
	  t += dt;
	  New_ket=ketstate;

	}
    }
  
  // double time=0.0;
  // for(int i = 0;i<N;i++)
  // {
  
  //   std::cout<<time/Tmax<<"  "<<X_avg[i]<<"  "<<Y_avg[i]<<"  "<<Z_avg[i]<<std::endl;
  //   time +=dt;
  // }
  return 0;
}
// --------------- Implementing functions ----------------
double Lambda(double time)
{
  return G*(1.0/cosh(Nu*(1.0 - (time/Tau))));
}

std::complex<double> Correlated_noise (double Noise1, double Noise2)
{
  arma::vec Eigenvalues;
  arma::mat Eigenvectors;
   // u that generates the noise 
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
  tweaked= (Q*original + target_mean);
  std::complex<double> Cor_Noise (tweaked(0),tweaked(1));
  return Cor_Noise;
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
arma::cx_mat Euler(arma::cx_mat & VectorState,double time,double dTime,std::complex<double> dW)
{
  arma::cx_mat Id(2,2);
  Id.eye();
  std::complex<double> Iimag(0.0,1.0);
  arma::cx_mat A = (Id - Jdagga*J*0.5 - (Iimag/Hbar)*Hamiltonian(time))*VectorState;
  arma::cx_mat B = J*VectorState;
  arma::cx_mat eulermethod = A*dTime+B*dW;
  VectorState += eulermethod;
  // return eulermethod;
  return VectorState;
}
arma::cx_mat State(arma::cx_mat VectorState)
{
  arma::cx_mat rho = VectorState*VectorState.t();
  return rho;
}
arma::cx_mat Normalise(arma::cx_mat Not_normalised)
{
  std::complex<double> Norm=std::sqrt(arma::as_scalar(Not_normalised.t()*Not_normalised));
  arma::cx_mat normalisedstate=Not_normalised/Norm;
  return normalisedstate;
}
arma::cx_mat Unitary(arma::cx_mat & VectorState,double time,double dTime, std::complex<double> dW)
{
  std::complex<double> Iimag(0.0,1.0);
  arma::cx_mat A = -(Iimag/Hbar)*Hamiltonian(time)*VectorState;
  VectorState += A*dTime;
  return VectorState;
}
arma::cx_mat Non_Unitary(arma::cx_mat & VectorState,double time,double dTime,std::complex<double> dW)
{
  arma::cx_mat Id(2,2);
  Id.eye();
  arma::cx_mat A = (Id - Jdagga*J*0.5)*VectorState;
  arma::cx_mat B = J*VectorState;
  VectorState +=  A*dTime + B*dW;
  return VectorState;
}
arma::cx_mat Initial_Conditions(std::complex<double> x1,std::complex<double> x2)
{
  arma::cx_mat Aux=Ket(x1,x2);
  return Aux;
}
