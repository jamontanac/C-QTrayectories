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

- [X] Implement the noise with normal distribution with parameters given by ![equation](https://latex.codecogs.com/gif.latex?\mu=0) and ![equation](https://latex.codecogs.com/gif.latex?\sigma=\sqrt{dt}). [Generate Normal noise](https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform)

- [X] Make use of the Mersenne twister engine by Matsumoto and Nishimura to generate the numbers, it's useful due to its period. the way to use it is for the case of a uniform distribution
	* the way to use it is for the case of a uniform distribution is
	```c++
	seed=std::random_device{}(); //Here we define the time as the seed
	std::mt19937 generator(seed);
	std::uniform_int_distribution<int> distribution(-100,100);// this produce a uniform distribution between -100 and 100
	```
	* And for the case of a normal distribution
	```c++
	seed=1;
	Mu=0.0;
	Sigma=1.0;
	std::mt19937 generator(seed);
  	std::normal_distribution<double> Noise(Mu, Sigma);//set a normal distribution to call it we do Noise(generator)
	```
	* if we want to improve the generator to 64 bit all we need is to call it as "mt19937_64"
	```c++
	seed=std::random_device{}();
	std::mt19937_64 generator(seed);
	```

- [X] Implement the normalisation function to normalise at every step of the method.

	* Up to this point everything is implemented in the code called:
	``` c++
	Secondpart.cpp
	```
	And the way I am compiling until now is:
	``` c++
	g++-8 -std=c++11 Secondpart.cpp -fsanitize=address -fsanitize=leak -fvisibility=hidden -larmadillo
	```

- [ ] Implement the code for a single trajectory in the range of ![equation](https://latex.codecogs.com/gif.latex?N=t_{max}/dt) with ![equation](https://latex.codecogs.com/gif.latex?t_{max}) the time of the experiment.

- [ ] At the level of a single trajectory compute the coordinates ![equation](https://latex.codecogs.com/gif.latex?x(t),y(t),z(t)) of the block sphere as ![equation](https://latex.codecogs.com/gif.latex?\mathrm{tr}(\sigma_{i}\rho(t))).

- [X] For a single trajectory compute the thermodynamic quantities such as heat ![equation](https://latex.codecogs.com/gif.latex?\mathcal{Q}[\rho(t)]) and work ![equation](https://latex.codecogs.com/gif.latex?\mathcal{W}[\rho(t)]) at every time.

- [ ] Change the method to a method that makes the Unitary and non-unitary evolution of a state in order to calculate separate things like heat, work and energy.

- [ ] Generate correlated noise to change the unravelings.[General Theory](https://es.wikipedia.org/wiki/Distribuci%C3%B3n_normal_multivariante)
> To generate a multivariable normal distribution it is possible to take a random vector of independent variables ![equation](https://latex.codecogs.com/gif.latex?X_1,X_2\dots,X_n) with normal distribution ![equation](https://latex.codecogs.com/gif.latex?X_{j}\stackrel{d}{=}\mathcal{N}(\mu_{j},\sigma_{j}^{2})) for ![equation](https://latex.codecogs.com/gif.latex?j=1,\dots,n) and then consider the next linear combination ![equation](https://latex.codecogs.com/gif.latex?Y=\sum_{j=1}^{n}\alpha_{j}X_j) we get that the random variable ![equation](https://latex.codecogs.com/gif.latex?Y) has a multivariate normal distribution with parameters ![equation](https://latex.codecogs.com/gif.latex?\mu:=\sum_{j=1}^n\mu_j\alpha_j) and ![equation](https://latex.codecogs.com/gif.latex?\sigma^2:=\sum_{j=1}^n\sigma_j^2\alpha_j^2). This for the case of independent variables the general case is presented in the section 5 o this [Book](http://bdigital.unal.edu.co/48054/2/9587014499.PDF).




