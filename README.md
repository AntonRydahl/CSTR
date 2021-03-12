Parallel Implementation of the CSTR Model
=========================================
This example demonstrates how to simulate the model dynamics of the Continuous Stirred Tank Reactor model for an example input profile.
This code implements the model given in the paper <a href="https://orbit.dtu.dk/en/publications/nonlinear-model-predictive-control-for-an-exothermic-reaction-in-">Nonlinear Model Predictive Control for an Exothermic Reaction in an Adiabatic CSTR</a>.

The simulations in this example are run in parallel using the Open Multi-processin (OpenMP).

Getting Started
---------------
Start out by navigating to <https://thinlinc.gbar.dtu.dk/main/>. In the top left corner, click *Applications* and select *DTU*. In order to get the result nicely displayed, open *xterm (VirtualGL-application-node)*. If you choose to use another terminal on thinlinc, the Matlab driver will not work and the results will not be displayed.

Create a new directory where you can run the code. For instance, type
```
mkdir 02686_CSTR
```
and navigate to the new directory
```
cd 02686_CSTR
```
In order to download these files, type
```
git init
git clone https://github.com/AntonRydahl/CSTR
```
and navigate to the downloaded files
```
cd cstr
```
You are now ready to run the example. In your folder you now have libraries containing a random number generator, *Mersenne Twister*, an implicit first order ODE solver, *Implicit Euler* and a library for generating a scalar Standard Wiener Process. The Newton solver uses in the Implicit Euler method uses LAPACK. You can eventually take a look at the makefile to see how the Fortran version of LAPACK can be used. 

The folder also contains a Matlab driver to illustrate the results.

Commands to Run the Driver
---------------------------
The driver can be executed with the command:
```
./driver.sh <number of realisations> <number of threads>
```

Expected Result
---------------
The C function *project.c* will compute a user specified number of realisations of model noise and simulate the model for the different realisations of noise. Afterwards the Matlab driver *driver.m* illustrates the solution.

![alt text](https://github.com/AntonRydahl/CSTR/blob/main/implicit_explicit.png) 
