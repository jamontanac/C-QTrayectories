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
5. The output of this code is redirect to the folder 