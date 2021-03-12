#include "ImplicitEulerSolver.h"
#include <stdbool.h>
#include <math.h>
#include <stdio.h>

void vector_implicit_euler(
    int N,
    int n,
    int NS,
    double *pt,
    double *px,
    double *pdW,
    double *workspace_lf,
    int *workspace_d,
    int max_iterations,
    double tolerance,
    functiontype f_func,
    functiontype g_func,
    functiontype J_func,
    double *pu,
    double *pd,
    void *pP,
    double *px0
){
    unsigned int row, col, sim, i,j,k;
    double h; // temporal step
    double *pF = &workspace_lf[0];
    double *pG = &workspace_lf[n];
    double *ppsi = &workspace_lf[n+n];
    double *workspace_inner = &workspace_lf[n+n+n];
    int len2dims_X = (N+1)*n;
    int len2dims_dW = len2dims_X-n;
    int index_X = 0;
    int index_dW = 0;

    for (sim = 0; sim < NS;sim++){
        // Imposing initial condition
        i = index_X;
        for (row = 0; row < n; row++) {
            px[i] = px0[row];
            i++;
        }

        // Initializing indexes
        i = index_X;
        j = index_X+n;
        k = index_dW;
        for (col = 0; col < N; col++) {
            // Invoking drift and diffusion terms
            f_func(&pt[col],&px[i],pu,pd,pP,pF);
            g_func(&pt[col],&px[i],pu,pd,pP,pG);

            // Calculating time step
            h = pt[col+1]-pt[col];

            // Initial guess for Newton solver
            for (row = 0; row < n; row++) {
                ppsi[row] = px[i] + (pG[row]*pdW[k]);
                px[j] = ppsi[row]+ (h*pF[row]);
                i++;
                j++;
                k++;
            }

            // Invoking newton solver
            newton_solver(
                f_func,
                J_func,
                max_iterations,
                tolerance,
                n,
                h,
                &pt[col],
                &px[i], // since i has been incremented to what was j
                ppsi,
                workspace_inner,
                workspace_d,
                pu,
                pd,
                pP
            );
        }
        index_X += len2dims_X;
        index_dW += len2dims_dW;
    }
}

void newton_solver(
    functiontype f_func,
    functiontype J_func,
    int max_iterations,
    double tolerance,
    int n,
    double dt,
    double *pt,
    double *px,
    double *ppsi,
    double *workspace_lf,
    int *workspace_d,
    double *pu,
    double *pd,
    void *pP
){
    double *pftemp = &workspace_lf[0];
    double *pjacobian = &workspace_lf[n];
    double *pdRdX = &workspace_lf[n*(n+1)];
    double *pR = &workspace_lf[n*(n+n+1)];

    f_func(pt, px,pu,pd,pP,pftemp);
    J_func(pt, px,pu,pd,pP,pjacobian);
    
    // Initializing residuals
    int i = 0;
    for (i=0;i<n;i++){
        pR[i] =  px[i]-pftemp[i]*dt-ppsi[i];
    }

    // Parameters for DGESV
    int N = n;
    int NRHS = 1;
    int LDA = N;
    int LDB = N;
    int INFO;
    // Minimizing residuals
    int iterations = 0; // of iterations

    // Boolean to check if the residuals meet the tolerance
    bool has_converged;

    // Index to check when a diagonal element is reached 
    int diagonal_index = 0;

    // Size of square matrix
    int n_square = n*n;
    for (iterations = 0;iterations<max_iterations;iterations++){
        // Initializing system matrix to solve
        for (i=0;i<n_square;i++){
            if (i == diagonal_index){ // Handling diagonal elements
                pdRdX[i] = 1 - pjacobian[i]*dt;
                diagonal_index += n+1;
            }
            else {
                pdRdX[i] = -pjacobian[i]*dt;
            }
        }

        // Minimizing residuals
        dgesv_(
            &N, 
            &NRHS,
            pdRdX,
            &LDA, 
            workspace_d,
            pR,
            &LDB,
            &INFO
        );

        // Updating the solution x
        for (i=0;i<n;i++){
            px[i] -= pR[i];
        }

        // Updating the residuals and calculating the infinity norm
        f_func(pt, px,pu,pd,pP,pftemp);
        has_converged = true;
        for (i=0;i<n;i++){
            pR[i] =  px[i]-pftemp[i]*dt-ppsi[i];
            has_converged &= fabs(pR[i]) < tolerance;
        }

        // If the infinity norm is less than the tolerance, the method terminates
        if (has_converged){
            break;
        }

        // Saving one call to the jacobian if convergence is obtained
        J_func(pt, px,pu,pd,pP,pjacobian);

        // Resetting the diagonal index
        diagonal_index = 0;
    }
}

void vector_implicit_euler_final_step(
    int N,
    int n,
    int NS,
    double *pt,
    double *px,
    double *pdW,
    double *workspace_lf,
    int *workspace_d,
    int max_iterations,
    double tolerance,
    functiontype f_func,
    functiontype g_func,
    functiontype J_func,
    double *pu,
    double *pd,
    void *pP,
    double *px0
){
    unsigned int row, col, sim;
    unsigned int k = 0;
    double h; // temporal step
    double *pF = &workspace_lf[0];
    double *pG = &workspace_lf[n];
    double *ppsi = &workspace_lf[n+n];
    double *pxtemp_1 = &workspace_lf[n+n+n]; // temporary x values for current step
    double *pxtemp_2 = &workspace_lf[n+n+n+n]; // temporary x values for last step
    double *workspace_inner = &workspace_lf[n+n+n+n+n];
    double *pxtemp; // used to switch current x step to previous x step

    int result_index = 0;
    for (sim = 0; sim < NS;sim++){
        // Running Euler Maruyama algorithm
        for (col = 0; col < N; col++) {
            // Invoking drift and diffusion terms
            if (col >0 ){
                f_func(&pt[col],pxtemp_1,pu,pd,pP,pF);
                g_func(&pt[col],pxtemp_1,pu,pd,pP,pG);
            }
            // Imposing initial condition
            else {
                f_func(&pt[col],px0,pu,pd,pP,pF);
                g_func(&pt[col],px0,pu,pd,pP,pG);
            }

            // Calculating time step
            h = pt[col+1]-pt[col];
            
            // There are now three cases: 
            // 1. Either we impose the initial condition stored in px0 
            if (col == 0){
                for (row = 0; row < n; row++) {
                    ppsi[row] = px0[row] + (pG[row]*pdW[k]);
                    pxtemp_2[row] = ppsi[row]+ (h*pF[row]);
                    k++;
                }

                // Invoking newton solver
                newton_solver(
                    f_func,
                    J_func,
                    max_iterations,
                    tolerance,
                    n,
                    h,
                    &pt[col],
                    pxtemp_2, // the result will be written here
                    ppsi,
                    workspace_inner,
                    workspace_d,
                    pu,
                    pd,
                    pP
                );

                // Updating pointers such that current time step is the previous time step
                pxtemp = pxtemp_2;
                pxtemp_2 = pxtemp_1;
                pxtemp_1 = pxtemp;
            }
            // 2. We iterate over all temporary solutions which are not stored
            else if (col < N-1){
                for (row = 0; row < n; row++) {
                    ppsi[row] = pxtemp_1[row] + (pG[row]*pdW[k]);
                    pxtemp_2[row] = ppsi[row]+ (h*pF[row]);
                    k++;
                }

                // Invoking newton solver
                newton_solver(
                    f_func,
                    J_func,
                    max_iterations,
                    tolerance,
                    n,
                    h,
                    &pt[col],
                    pxtemp_2, // the result will be written here
                    ppsi,
                    workspace_inner,
                    workspace_d,
                    pu,
                    pd,
                    pP
                );

                // Updating pointers such that current time step is the previous time step
                pxtemp = pxtemp_2;
                pxtemp_2 = pxtemp_1;
                pxtemp_1 = pxtemp;
            }
            // 3. Or we store the final solution computed in step N+1
            else {
                for (row = 0; row < n; row++) {
                    ppsi[row] = pxtemp_1[row] + (pG[row]*pdW[k]);
                    px[result_index] = ppsi[row]+ (h*pF[row]);
                    k++;
                    result_index++;
                }

                // Invoking newton solver
                newton_solver(
                    f_func,
                    J_func,
                    max_iterations,
                    tolerance,
                    n,
                    h,
                    &pt[col],
                    &px[result_index-n], // The result is written to the final solution
                    ppsi,
                    workspace_inner,
                    workspace_d,
                    pu,
                    pd,
                    pP
                );
            }
        }
    }
}
