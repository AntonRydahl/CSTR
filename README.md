Parallel Implementation of the CSTR Model
=========================================
This example demonstrates how to simulate the model dynamics of the Continuous Stirred Tank Reactor model for an example input profile. The simulations are run in parallel using the Open Multi-processin (OpenMP).

Getting Started
---------------
Start out by navigating to <https://thinlinc.gbar.dtu.dk/main/>. In the top left corner, click *Applications* and select *DTU*. In order to get the result nicely displayed, open *xterm (VirtualGL-application-node)*.

Create a new directory where you can run the code. For instance, type
´´´
mkdir 02686_CSTR
´´´
Navigate to the new directory:
´´´
git init
git clone https://github.com/AntonRydahl/CSTR
´´´
You are now ready to run the example. In your folder you now have libraries containing a random number generator, *Mersenne Twister*, an implicit first order ODE solver, *Implicit Euler* and a library for generating a scalar Standard Wiener Process.

The folder also contains a Makefile for the project and a Matlab driver to illustrate the results.

Commands to Run the Example
---------------------------
The driver can be executed with the command:
```
$ ./driver.sh <number of realisations> <number of threads>
```

Expected Result
---------------
The C function *parallel_driver.c* will compute a user specified number of realisations of model noise and simulate the model for the different realisations of noise. Afterwards the Matlab driver *driver.m* illustrates the solution.

![alt text](https://github.com/AntonRydahl/CSTR/implicit_explicit.png) 
