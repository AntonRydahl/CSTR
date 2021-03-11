Parallel Implementation of the CSTR Model
=========================================
This example demonstrates how to simulate the model dynamics of the CSTR model for an example input profile. The implementation will use the defualt number of logical CPUs on the system. Please make sure that OpenMP is installed before running the makefile.

Command to Run the Example
---------------------------
The driver can be executed with the command:
```
$ ./driver.sh <number of realisations>
```

Expected Result
---------------
The C function *parallel_driver.c* will compute a user specified number of realisations of model noise and simulate the model for the different realisations of noise. Afterwards the Matlab driver *driver.m* illustrates the solution.

![alt text](https://github.com/AntonRydahl/stochastic_differential_equations/blob/main/CSTR/implicit_parallel/implicit_explicit.png) 
