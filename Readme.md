# Quantum Trajectories in C++
In the present repository I present an implementation of a quantum trajectory in the case of diffusive measurements.
The library used here are:
* ### **Complex** [Info](http://www.cplusplus.com/reference/complex/)
* ### **Armadillo** [Info](http://arma.sourceforge.net/)
* ### **Openmp** [Info](https://www.openmp.org/)
And the way of compiling is simply g++ -std=c++11 "Code.cpp" -larmadillo
To configure the GitHub i used the help of [Don't be afraid to commit](http://dont-be-afraid-to-commit.readthedocs.io/en/latest/git/commandlinegit.html)
To do list:
- [X] Implement as constants the part of the operators such as the Pauli matrix operators (![equation](https://latex.codecogs.com/gif.latex?$\sigma_i$)) and the operators of the bath (![equation](https://latex.codecogs.com/gif.latex?J)) and (![equation](https://latex.codecogs.com/gif.latex?$J^{\dag}$))