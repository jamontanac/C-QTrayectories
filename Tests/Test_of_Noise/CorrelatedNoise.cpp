#include <iostream>
#include <cmath>
#include <complex>
// #include "armadillo"
#include "omp.h"
int main(void)
{
  std::mt19937 generator(1);
  std::normal_distribution<double> Noise1(Mu, Sigma);
  std::normal_distribution<double> Noise2(Mu, Sigma);
  return 0;
}
