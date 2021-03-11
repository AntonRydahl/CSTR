/**
* @snippet CSTR.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "CSTR.h"
#include <math.h>
#include "ImplicitEulerSolver.h"

void flow_rate(double *parray){
    parray[0] = 700;
    parray[1] = 700;
    parray[2] = 700;
    parray[3] = 600;
    parray[4] = 600;
    parray[5] = 500;
    parray[6] = 500;
    parray[7] = 400;
    parray[8] = 400;
    parray[9] = 300;
    parray[10] = 300;
    parray[11] = 300;
    parray[12] = 200;
    parray[13] = 200;
    parray[14] = 200;
    parray[15] = 200;
    parray[16] = 300;
    parray[17] = 300;
    parray[18] = 400;
    parray[19] = 400;
    parray[20] = 500;
    parray[21] = 500;
    parray[22] = 600;
    parray[23] = 600;
    parray[24] = 700;
    parray[25] = 700;
    parray[26] = 700;
    parray[27] = 700;
    parray[28] = 200;
    parray[29] = 200;
    parray[30] = 200;
    parray[31] = 200;
    parray[32] = 700;
    parray[33] = 700;
    parray[34] = 700;
}

CSTR_parameters default_parameters(){
    struct CSTR_parameters params;
    params.final_time = 35*60;
    params.EaR = 8500;
    params.rho = 1;
    params.DeltaH = -560000;
    params.cP = 4186;
    params.beta = 133.7793;
    params.CAin = 0.8;
    params.CBin = 1.2;
    params.Tin = 273.65;
    params.V = 0.105;
    params.k0 = 48266327438.6281;
    params.sigma = 10;
    return params;
}

void CSTR_3D_drift(
    double *pt,
    double *px, 
    double *pu, 
    double *pd, 
    void *pP,
    double *pxdot
){
    CSTR_parameters *params = (CSTR_parameters *) pP;
    double CA = px[0]; // Concentration of compound A 
    double CB = px[1]; // Concentration of compound B 
    double temperature = px[2]; // Temperature named T in the model
    double f = pu[0]; // scaling to [seconds / 10000]

    //  Arrhenius expression
    double k_arrhenius = params->k0*exp(params->EaR*(-1/temperature));
    double r = k_arrhenius*CA*CB;

    // Production rate and rate of change in temperature
    double RA = -r;
    double RB = -2*r;
    double RT = params->beta*r;

    double FV = f/params->V;

    // Derivatives
    pxdot[0] = FV*(params->CAin-CA) + RA; // CA_bar
    pxdot[1] = FV*(params->CBin-CB) + RB; //CB_bar
    pxdot[2] = FV*(params->Tin-temperature) + RT; // T_bar
}

void CSTR_3D_diffusion(
    double *pt,
    double *px, 
    double *pu, 
    double *pd, 
    void *pP,
    double *pxdot
){
    CSTR_parameters *params = (CSTR_parameters *) pP;
    double f = pu[0]; // scaling to [seconds / 10000]
    pxdot[0] = 0;
    pxdot[1] = 0;
    pxdot[2] = (params->sigma*f)/params->V;
} 

void CSTR_3D_drift_jacobian(
    double *pt,
    double *px, 
    double *pu, 
    double *pd, 
    void *pP,
    double *pxdot
){
    CSTR_parameters *params = (CSTR_parameters *) pP;
    double CA = px[0]; // Concentration of compound A 
    double CB = px[1]; // Concentration of compound B 
    double temperature = px[2]; // Temperature named T in the model
    double f = pu[0]; // scaling to [seconds / 10000]

    //  Arrhenius expression
    double k_arrhenius = params->k0*exp(params->EaR*(-1/temperature));

    double FV = f/params->V;
    double kCA = k_arrhenius*CA;
    double kCB = k_arrhenius*CB;
    double kT = CA*CB*params->EaR*params->k0*exp(params->EaR*(-1/temperature));
    kT /= (temperature*temperature);
    // It may look like the jacobian is implemented as the transpose, but it simply uses col major storage
    pxdot[0] = -FV-kCB;
    pxdot[1] = -(kCB+kCB);
    pxdot[2] = params->beta*kCB;

    pxdot[3] = -kCA;
    pxdot[4] = -FV-(kCA+kCA);
    pxdot[5] = params->beta*kCA;

    pxdot[6] = kT;
    pxdot[7] = -(kT+kT);
    pxdot[8] = -FV+params->beta*kT;
}


void implicit_simulation(
    double *pt,
    double *px,
    double *pdW,
    double *pworkspace_lf,
    int *pworkspace_d,
    int max_iterations,
    double tolerance,
    functiontype f_func,
    functiontype g_func,
    functiontype J_func,
    double *pu,  // for storing closed loop input profiles
    double *pd,
    CSTR_parameters *pP,
    int num_realizations,
    int num_samples,
    int time_steps_per_sample,
    int N,
    int n,
    int dw_increment, // if 0, the same measurement noise will be used in all simulations
    int p_increment // if 0, the same parameter vector will be used in all simulations
){
    int i = 0;
    int j = 0;
    int x_index_outer = 0;
    int x_index_inner = 0;
    int dw_index_outer = 0;
    int dw_index_inner = 0;
    int t_index = 0;
    int p_index = 0;
    int x_increment = n*(N+1);
    int sample_size = n*time_steps_per_sample;
    for (i=0;i<num_realizations;i++){
        x_index_inner = x_index_outer;
        dw_index_inner = dw_index_outer;
        for (j=0;j<num_samples;j++){
            vector_implicit_euler(
                time_steps_per_sample,
                n,
                1,
                &pt[t_index],
                &px[x_index_inner],
                &pdW[dw_index_inner],
                pworkspace_lf,
                pworkspace_d,
                max_iterations,
                tolerance,
                f_func,
                g_func,
                J_func,
                &pu[j],
                pd,
                &pP[p_index],
                &px[x_index_inner]
            );
            x_index_inner += sample_size;
            dw_index_inner += sample_size;
            t_index += time_steps_per_sample;
        }
        x_index_outer += x_increment;
        dw_index_outer += dw_increment;
        t_index = 0;
        p_index += p_increment;
    }
}
