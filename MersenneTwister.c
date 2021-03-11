/// @file MersenneTwister.c

#include "MersenneTwister.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static unsigned long mersenne_seed = 4357;

void set_seed(
    unsigned int x
){
  mersenne_seed = x;
};

/*******************************************************************************
Helper functions for Mersenne Twister
*******************************************************************************/

// Initializing the generator
void seed_mt(
    unsigned long *parr
) {
    // Length of generator
    int n = 624;
    // 0xffffffff is equivalent to 4294967295 in decimal
    // 0xffffffffUL & is the modulo operaion in the field with 2^32 elements
    parr[0] = 0xffffffffUL & mersenne_seed;
    int i = 0;
    for (i=1;i<n;i++) {
        parr[i] = 0xffffffffUL & (69069 * parr[i-1]);
    }
}

// Twisting the generator
void twist(
  unsigned long *parr
){
    // Length of generator
    int n = 624;
    // Offset between indexes in generator
    int m = 397;
    // Iterators
    int i = 0;
    int j,k;
    // Random integer
    unsigned int x;
    for (i=0;i<n;i++){
        j = (i+1);
        if (j>=n){
            j=0;
        }
        k = (i + m);
        if (k>=n){
            k=0;
        }
        // Concatenation of lower bits and upper bits
        x = (parr[i] & 0x80000000UL) | (parr[j] & 0x7fffffffUL);
        if (x & 0x1){ // If x is odd, xor with A, which is 0x9908b0dfUL
            // Fast multiplication with Frobenius matrix
            parr[i] = parr[k] ^ (x >> 1) ^ 0x9908b0dfUL;
        }
            else { // If x is even, don't
            parr[i] = parr[k] ^ (x >> 1);
        }
    }
}

/*******************************************************************************
Implementation of Mersenne Twister
*******************************************************************************/

void mersenne_twister(
    double *parr, 
    unsigned long *pgenerator, 
    int N
){
    int i;
    // Length of generator
    int n = 624;
    // Offset between indexes in generator
    seed_mt(pgenerator);
    twist(pgenerator);
    int index = 0;
    for(i=0;i<N;i++){
        // If all values in generator has been used, update the generator
        if (index >= n) {
            twist(pgenerator);
            index = 0;
        }
        // Tempering operation
        unsigned long y = pgenerator[index];
        y = y ^ (y >> 11);
        y = y ^ ((y << 7) & 0x9d2c5680);
        y = y ^ ((y << 15) & 0xefc60000);
        y = y ^ (y >> 18);
        index += 1;
        // Transforming unsigned integer back to the open unit interval (0,1)
        parr[i] = (0.5+(long double) y )/(0xffffffffUL+1);
    }
    mersenne_seed = pgenerator[n-1];
}

/*******************************************************************************
Box Muller Transformation of Uniform Distribution
*******************************************************************************/

void box_muller(
    double *parr,
    int N,
    long double mu, 
    long double sigma
){
    int i;
    long double U0,U1,Z0,Z1,R,Omega,firstVal;
    long double PI = 3.14159265;
    // If the length is uneven, the first value is
    // also used as the last
    firstVal = parr[0];
    for(i=0;i<N;i++){
        U0 = parr[i];
        U1 = (i < N-1) ? parr[i+1] : firstVal;
        R = sqrt(-2*log(U0));
        Omega = 2*PI*U1;
        Z0 = R*cos(Omega);
        Z1 = R*sin(Omega);
        parr[i] = mu+sigma*Z0;
        if (i<N-1){
            parr[i+1] = mu+sigma*Z1;
            i++;
        }
    }
}

void d_rand_standard_normal(
    double *parr, 
    unsigned long *pworkspace, 
    int N
){
  mersenne_twister(parr,pworkspace,N);
  box_muller(parr,N,0.0,1.0);
}

void d_rand_normal(
    double *parr, 
    unsigned long *pworkspace, 
    int N,
    long double mu, 
    long double sigma
){
  mersenne_twister(parr,pworkspace,N);
  box_muller(parr,N,mu,sigma);
}
