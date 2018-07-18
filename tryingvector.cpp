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
    std::cout<<A[i]<<"  ";
  }
  // std::cout<< A<<std::endl;
  return 0;

}
