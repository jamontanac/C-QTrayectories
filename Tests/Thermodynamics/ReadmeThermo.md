## Thermodynamics Readme
In the present folder i present the solution of the first part of computing the quantum thermodynamic quantities.

* The code named [NewImplementation](./NewImplementation.cpp) contains the implementation of the code that computes:

1. The components of the Bloch Sphere ![equation](https://latex.codecogs.com/gif.latex?x,&space;y,&space;z) in function of the time.

2. The Work, Heat and Energy of the system by evolving unitarily and non unitarily and saving these evolutions to compute these quantities.

3. It also compute the averages over trajectories and save it in a vector. And the way of printing these quantities is as follows
	``` c++
	double time=0.0;
    	for(int i = 0;i<N;i++)
    	{
    		std::cout<<time/Tmax<<"  "<<X_avg[i]<<"  "<<Y_avg[i]<<"  "<<Z_avg[i]<<std::endl;
    		time +=dt;
    	}
	```
4. The way I am compiling in order to check any kind of errors or leaks in the memory or maybe variables that are declared but are not being used, and the way to do this is as follows:
	```c++
	g++-8 -std=c++11 NewImplementation.cpp -larmadillo -fsanitize=address -fsanitize=leak -fsanitize=undefined -fvisibility=hidden -Wall
	```
5. The output of this code is redirect to the folder [Plots](./Plots) in three different files named:

	* **datos.txt**
	* **datosAVG.txt**
	* **thermo.txt**

In the first one is the data of the Bloch sphere over a single trajectory. The second one is the The Bloch sphere for averaged trajectories (The number of trajectories that I use to average can change). Finally the third one is the data of the thermodynamic quantities along a single trajectory.

* The folder named [Plots](./Plots) has the files to plot de output given by Newimplementation.cpp, there are some files named:
	* **plotsAVG.py**
	* **Thermoplots.py**
	* **Blochplots.py**
These three are the programs in python used to plot the averaged trajectories, the thermodynamic quantities and the Bloch sphere respectly. These generates the .pdf files that can be found in the folder [Plots](./Plots).

* The File named [Newimplementation_CorNoise](./Newimplementation_CorNoise.cpp) is an implementation of the last code but with the difference that this one use correlated noise based in the noise produced by a normal distribution (![equation](https://latex.codecogs.com/gif.latex?\mu = 0, &space; \sigma = 1)). 
	- In the program to change the unravelings we only have to change the ![equation](https://latex.codecogs.com/gif.latex?\U_{corr})
		```c++
			// u that generates the noise 
  			std::complex<double> U_corr (1.0,0.0);
		```
		
	By changing this ![equation](https://latex.codecogs.com/gif.latex?\U_{corr}) the noise I generate changes for example:
	> ![equation](https://latex.codecogs.com/gif.latex?U_{corr} =1), generates real noise (Homodyne Detection), ![equation](https://latex.codecogs.com/gif.latex?U_{corr}=-1), generates imaginary pure noise (Heterodyne Detection) and ![equation](https://latex.codecogs.com/gif.latex?U_{corr}=0) generates complex non correlated noise (Photo detection)
		