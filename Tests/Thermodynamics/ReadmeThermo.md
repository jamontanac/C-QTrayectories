## Thermodynamics Readme
In the present folder i present the solution of the first part of computing the quantum thermodynamic quantities.

* The code named [NewImplementation](./NewImplementation.cpp) contains the implementation of the code that computes:

1. The components of the Bloch Sphere ![equation](https://latex.codecogs.com/gif.latex?x,&space;y,&space;z) in function of the time.

2. The Work, Heat and Energy of the system by evolving unitarily and non unitarily and saving these evolutions to compute these quantities.

3. It also compute the averages over trajectories and save it in a vector.
	``` c++
	double time=0.0;
    	for(int i = 0;i<N;i++)
    	{
    		std::cout<<time/Tmax<<"  "<<X_avg[i]<<"  "<<Y_avg[i]<<"  "<<Z_avg[i]<<std::endl;
    		time +=dt;
    	}
	```
