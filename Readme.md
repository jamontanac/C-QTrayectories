# Quantum Trajectories in C++
In the present repository I present an implementation of a quantum trajectory in the case of diffusive measurements.
The library used here are:
* ### **Complex** [Info](http://www.cplusplus.com/reference/complex/)
* ### **Armadillo** [Info](http://arma.sourceforge.net/)
* ### **Openmp** [Info](https://www.openmp.org/)
And the way of compiling is simply g++ -std=c++11 "Code.cpp" -larmadillo
To configure the GitHub i used the help of [Don't be afraid to commit](http://dont-be-afraid-to-commit.readthedocs.io/en/latest/git/commandlinegit.html)
To do list:

- [X] Implement as constants the part of the operators such as the Pauli matrix operators (![equation](https://latex.codecogs.com/gif.latex?$\sigma_i$)) and the operators of the bath (![equation](https://latex.codecogs.com/gif.latex?J)) and (![equation](https://latex.codecogs.com/gif.latex?J^{\dag})).

- [X] Implement the functions of the problem such as the driving ![equation](https://latex.codecogs.com/gif.latex?\lambda(t)), the Hamiltonian, the Ket and Bra operators as well as the State function which returns the product of Ket and Bra to generate the density matrix ![equation](https://latex.codecogs.com/gif.latex?\rho(t)).

- [X] Implement the Stochastic  Euler method to solve the stochastic differential equation. 

- [X] Implement the noise with normal distribution with parameters given by ![equation](https://latex.codecogs.com/gif.latex?\mu=0) and ![equation](https://latex.codecogs.com/gif.latex?\sigma=\sqrt{dt}).

- [X] Implement the normalisation function to normalise at every step of the method.

- [ ] Implement the code for a single trajectory in the range of ![equation](https://latex.codecogs.com/gif.latex?N=t_{max}/dt) with ![equation](https://latex.codecogs.com/gif.latex?t_{max}) the time of the experiment.
