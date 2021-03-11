#ifndef CSTR_MODEL_PARAMETERS
#define CSTR_MODEL_PARAMETERS

/**
 * This generec function type will be used throughout all solvers in this library. It is meant for returning \f$f\f$ in 
 * \f$dx = f(t,x,u,d,p) dt\f$ where \f$t\f$ is the temporal solution, \f$x\f$ is the spatial solution, \f$u\f$ is the control parameter,
 * \f$d\f$ is the disturbance and \f$p\f$ is an arbitrary parameter input. 
 * The function can be impemented to return a vector or a matrix dependent on the mathematical model and the solver. 
 * In case of a matrix, ensure you use column major storage of the matrix elements.
 *
 * The main purpose it is to allow for any kind of custom parameters in the function
 * calls. The functions of this type are meant for returning a spatial detivative 
 * which must be writeen to result.
 * @param[in] pt: Pointer to the temporal solution value(s).
 * @param[in] px: Pointer to the spatial solution value(s).
 * @param[in] pu: Pointer to the control parameter.
 * @param[in] pd: Pointer to the disturbances.
 * @param[in] pP: Void pointer to an arbitrary datastructure containing the custom parameters used in the function.
 * @param[in,out] pxdot: Pointer to the drift output \f$f\f$ (or equivalenty the diffusion).
 * 
 * @author Anton Rydahl
 * @date 22nd of October 2020
 */

typedef void (*functiontype)(
    double *pt,
    double *px, 
    double *pu, 
    double *pd, 
    void *pP, 
    double *pxdot
);

/**
 * The follwong function returns the flow rate F for a 35 minutes simulation.
 * The unit is [milliliter / minute]
 *
 * @param[in] parray: Pointer to array of doubles of length 35.
 *
 * @author: Anton Rydahl
 * @date: 8th of October 2020
 */

void flow_rate(double *parray);

/**
 * The follwong struct type contains all constants needed in the CSTR model.
 *
 *@author: Anton Rydahl
 *@date: 8th of October 2020
 */

typedef struct CSTR_parameters{
   double *flow_rate;
   double final_time;
   double EaR;
   double rho;
   double DeltaH;
   double cP;
   double beta;
   double CAin;
   double CBin;
   double Tin;
   double V;
   double k0;
   double sigma;
} CSTR_parameters;

/**
 * The following values and units are used for the parameters:
 * flow_rate:           [mL / min]
 * final time: 35       [min]
 * EaR: 8500            [K]
 * rho : 1.0            [kg / L]
 * DeltaH: -560         [kJ / mol]
 * cP: 4.186           [(kJ)/(kg * K)]
 * beta: 133.7793
 * CAin: 0.8           [mol / L]
 * CBin: 1.2           [mol / L]
 * Tin: 273.65         [K]
 * V: 0.105             [L]
 * k0: 4.8266*10^10     [L / (mol * S)]
 * sigma: 10 
 *
 *@author: Anton Rydahl
 *@date: 8th of October 2020        
 */

CSTR_parameters default_parameters(); 

/**
 * Drift term for the 3 dimensional CSTR model. 
 * 
 * @param[in] pt: The temporal solution. Must be 1*sizeof(double).
 * @param[in] px: The spatial solution. Must contain \f$C_A\f$, \f$C_B\f$ and \f$T\f$ in this order. Must be 3*sizeof(double).
 * @param[in] pu: Pointer to the control parameter \f$F(t)\f$.
 * @param[in] pd: Pointer to disturbance. the input is unused in this case and can be set to NULL.
 * @param[in] pP: Void pointer to CSTR_parameters pointer containing all the constants used in the function.
 * @param[in,out] pxdot: The derivative of the function. Will contain the derivatives of \f$C_A\f$,
 * \f$C_B\f$ and \f$T\f$ in this order.
 *
 * @author: Anton Rydahl
 * @date: 8th of October 2020
 */

void CSTR_3D_drift(double *pt,double *px, double *pu, double *pd, void *pP,double *pxdot);

/**
 * Diffusion term for the 3 dimensional CSTR model. This function return only a three dimensional vector and hence
 * it is designed for the solver vectorEulerMaruyama().
 * 
 * @param[in] pt: The temporal solution. Must be 1*sizeof(double).
 * @param[in] px: The spatial solution. Must contain \f$C_A\f$, \f$C_B\f$ and \f$T\f$ in this order. Must be 3*sizeof(double).
 * @param[in] pu: Pointer to the control parameter \f$F(t)\f$.
 * @param[in] pd: Pointer to disturbance. the input is unused in this case and can be set to NULL.
 * @param[in] pP: Void pointer to CSTR_parameters pointer containing all the constants used in the function.
 * @param[in,out] pxdot: Will contain the disturbance of \f$C_A\f$,
 * \f$C_B\f$ and \f$T\f$ in this order.
 *
 * @author: Anton Rydahl
 * @date: 8th of October 2020
 */

void CSTR_3D_diffusion(double *pt,double *px, double *pu, double *pd, void *pP,double *pxdot);

/**
 * This function writes the Jacobian of the diffusion term to pxdot. The jacobain is meant for use in the solver vectorImplicitEuler().
 * The Jacobian is stored in column major order since the aforementioned solver used the Fortran routine DGESV.
 * 
 * @param[in] pt: The temporal solution. Must be 1*sizeof(double).
 * @param[in] px: The spatial solution. Must contain \f$C_A\f$, \f$C_B\f$ and \f$T\f$ in this order. Must be 3*sizeof(double).
 * @param[in] pu: Pointer to the control parameter \f$F(t)\f$.
 * @param[in] pd: Pointer to disturbance. the input is unused in this case and can be set to NULL.
 * @param[in] pP: Void pointer to CSTR_parameters pointer containing all the constants used in the function.
 * @param[in,out] pxdot: Will contain the Jacobain after execution. Must be of size 9*sizeof(double)-
 * \f$C_B\f$ and \f$T\f$ in this order.
 *
 * @author: Anton Rydahl
 * @date: 19th of November 2020
 */

void CSTR_3D_drift_jacobian(double *pt,double *px, double *pu, double *pd, void *pP,double *pxdot);



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
);


#endif
