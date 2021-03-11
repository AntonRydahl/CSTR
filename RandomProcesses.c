#include "MersenneTwister.h"
#include "RandomProcesses.h"
#include <stdio.h>

void scalar_wiener_process(
    double *pdW, 
    unsigned long *pworkspace,
    double FT,
    int N, 
    int nw,
    int NS
){
    // Defining time step
    double dt = (double) FT/(double)N;
    double sqrtdt = sqrt(dt);

    // Generating random numbers with mean 0 and variance dt
    d_rand_normal(pdW,pworkspace,nw*N*NS,0,sqrtdt);
}

void linspace(
    double *pT,
    double t0,
    double FT,
    int N
){
    int i = 0;
    double dt = (double) FT/(double)N;
    // Generating equidistant temporal vector
    for (i=0;i<(N+1);i++){
        pT[i] = t0;
        t0+= dt;
    }
}

void cumsum(
    double *pW,
    double *pdW,
    int N,
    int nw,
    int NS
){
    // Iterating dW
    int i = 0;
    // Iterating pW
    int j = 0;

    int total_size = N*nw*NS;
    int w_increment = (1+N)*nw;
    int next_simulation = 0;

    for (i=0;i<total_size;i++){
        if (j==next_simulation){
            for (j=j;j<next_simulation+nw;j++){
                pW[j] = 0;
            }
            next_simulation += w_increment;
        }
        pW[j] = pdW[i] + pW[j-nw];
        j++;
    }
}


