/**
* @snippet CSTR_driver.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "ImplicitEulerSolver.h"
#include "MersenneTwister.h"
#include "RandomProcesses.h"
#include "CSTR.h"

int main(int argc, char *argv[]){
    if (argc!=2){
        printf("Please provide the number of realizations of noise.\n");
        return 0;
    }


    // One time step is 10 seconds
    int time_steps_per_sample = 60;
    
    // Sample time is one minute
    int sample_time_seconds = 60;
    
    // Experiment takes 35 minutes
    int number_of_samples = 35;
    
    int total_steps = number_of_samples*time_steps_per_sample;

    // Dimensions used in solver
    // Number of time steps
    int N = total_steps;
    
    // Number of states in problem, concentration A, concentration B, temperature T
    int n = 3;
    
    // The dimension of the process noise is the same as the number of states
    int nw = 3;
    
    // Number of simulations of the experiment
    int NS = atoi(argv[1]); // Number of realizations of noise
    if (NS < 1){
        printf("Error: The number of simulations must be larger than 0.\n");
    }
    
    // Pointing to the drift
    functiontype f_func = CSTR_3D_drift;
    
    // Pointing to the diffusion
    functiontype g_func = CSTR_3D_diffusion;
    
    // Pointing to the Jacobian
    functiontype J_func = CSTR_3D_drift_jacobian;    

    // Allocating memory for the spatial solution
    int problem_size_x = n*(N+1)*NS;
    
    // The noise is shorter since there is no noise on the initial condition
    int problem_size_dW = nw*N*NS;
    double *pX = (double*) malloc(problem_size_x*sizeof(double));

    // Allocating memory for the white noise
    double *pdW = (double*) malloc(problem_size_dW*sizeof(double));

    // Allocating memory for the temporal solution
    double *pT = (double*) malloc((N+1)*sizeof(double));

    // Allocating memory for seeding Mersenne Twister
    unsigned long *pworkspace_ul = (unsigned long*) malloc(624*sizeof(unsigned long));

    // Allocating memory for the implicit-explicit Euler scheme for vector drift
    int max_num_threads = 8;
    double *pworkspace_lf = (double*) malloc(max_num_threads*(5+2*n)*n*sizeof(double));

    // Allocating memory for DGESV
    int workspace_d[max_num_threads*n];
    int *pworkspace_d = &workspace_d[0];

    // Allocating memory for the flow rate
    double *pflow_rate = (double*) malloc(number_of_samples*sizeof(double));

    ///! [allocating memory]

    ///! [parameters]
    // Declaring default parameters
    CSTR_parameters params= default_parameters();

    // Creating pointer to parameters
    CSTR_parameters *pP = &params;
    pP->sigma = 10;

    // There are no disturbances in this model
    double *pd = NULL;

    // Inserting default flow rate parameters for simulation
    flow_rate(pflow_rate);
    int i = 0;
    FILE *F_file = fopen("F.txt", "w");
    for (i=0;i<number_of_samples;i++){
        fprintf(F_file,"%1.15f\n",pflow_rate[i]);
        pflow_rate[i] = pflow_rate[i]/(60*1000);
    }
    fclose(F_file);

    // Imposing initial condition
    int x0_index = 0;
    int x0_increment = n*(N+1);
    for (i=0;i<NS;i++){
        pX[x0_index+0] = 0.05;
        pX[x0_index+1] = 0.25;
        pX[x0_index+2] = params.Tin;
        x0_index += x0_increment;
    }
    
    // Seeding the algorithm
    set_seed(12345);
    
    // Generating the noise - despite the name of the function it is noise.
    // It is not accumulated into a Brownian path
    scalar_wiener_process(
        pdW,
        pworkspace_ul,
        number_of_samples*sample_time_seconds, // to seconds
        N,
        nw,
        NS
    );
    
    // Generating equidistant time grid
    linspace(
        pT,
        0,
        number_of_samples*sample_time_seconds,
        N
    );

    // Parameters for implicit Euler solver
    int max_iterations = 20;
    
    // Tolerance for the Newton solver
    double tolerance = 10e-6;

    // The same parameters are used in every simulation
    // For monte Carlo simulations, set p_increment to 1
    // and generate a vector of parameters structs
    int p_increment = 0;

    // A new realization of noise is used in every simulation
    int dw_increment = nw*N;

    // For OpenMP loop
    int size_x = n*(N+1);

    // OpenMP thread private variables:
    int thread_index, number_of_threads, thread_points, thread_start;
    
    // Starting timing
    double timer = omp_get_wtime();
    
    #pragma omp parallel default(shared) private(thread_index,\
    number_of_threads, thread_points, thread_start)
    {   
        thread_index = omp_get_thread_num();
        number_of_threads =  omp_get_num_threads();
        thread_points = (int) NS / number_of_threads;
        thread_start = thread_index*thread_points;
        if (thread_index == number_of_threads-1){
            thread_points = NS - thread_start;
        }
        implicit_simulation(
            pT,
            &pX[thread_start*size_x],
            &pdW[thread_start*dw_increment],
            &pworkspace_lf[(2*n+5)*n*thread_index],
            &pworkspace_d[n*thread_index],
            max_iterations,
            tolerance,
            f_func,
            g_func,
            J_func,
            pflow_rate, //pu,  // for storing closed loop input profiles
            pd,
            pP,
            thread_points,
            number_of_samples,
            time_steps_per_sample,
            N,
            n,
            dw_increment, // if 0, the same measurement noise will be used in all simulations
            p_increment // if 0, the same parameter vector will be used in all simulations
        );
    }
    
    // Finishing timing
    timer = omp_get_wtime()-timer;
    printf("%lf\n", timer);
    
    
    FILE* X_file;
    FILE* T_file;
    X_file = fopen("X.txt", "w");
    T_file = fopen("T.txt", "w");
    int j,l;
    for (j=0;j<(N+1)*n*NS;j++){
        fprintf(X_file,"%1.15f\n",pX[j]);
    }
    for (l=0;l<(N+1);l++){
        fprintf(T_file,"%1.15f\n",pT[l]/60);
    }
    
    // Closing files
    fclose(T_file);
    fclose(X_file);
    
    // Avoiding memory leakage
    free(pflow_rate);
    free(pworkspace_lf);
    free(pworkspace_ul);
    free(pT);
    free(pdW);
    free(pX);

    return 0;
}

