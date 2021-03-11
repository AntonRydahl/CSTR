/// @file ImplicitEulerSolver.h

#ifndef EXPLICIT_IMPLICIT_EULER
#define EXPLICIT_IMPLICIT_EULER


///@cond IGNORE_IN_DOCUMENTATION
extern void dgemv_(char *trans, int *m, int *n,
              double *alpha, double *A, int *lda,
              double *px, int *incx, double *beta,
              double *Y, int *incy);

extern void dgesv_(int *N, int *NRHS,double *A,
              int *LDA, int *IPIV,double *B,
              int *LDB, int *INFO);
              
///@endcond

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
 * Implementation of the implicit-explicit Euler method. This numerical scheme approximates 
 * \f$dx(t) = f\big(x(t)\big)dt+g\big(x(t)\big)d\omega(t)\f$ by 
 * \f$x_{n+1} =x_{n} f(t_{n},x_{n})(t_{n+1}-t_{n})+g(t_{n},x_{n})d\omega_{n}\f$.
 * The implicit (backward) step in the diffusion term is solved with the Newton solver
 * newton_solver(). The jacobian of the diffusion term must be computed and stored in column
 * major order in J_func.
 * 
 * @param[in] N: The number of time steps (excluding the initial condition).
 * @param[in] n: The dimension of \f$x(t)\f$.
 * @param[in] NS: The number of simulations.
 * @param[in] pt: Pointer to the temporal solution.  Must be of size \f$(N+1)\cdot\text{sizeof}(\text{double})\f$.
 * @param[out] px: After operation this array will contain the spatial solution. Must be size \f$n\cdot (N+1)\cdot NS\cdot\text{sizeof}(\text{double})\f$.
 * @param[in] pdW: White noise. Must be size \f$n\cdot (N+1)\cdot NS \cdot\text{sizeof}(\text{double})\f$.
 * @param[in] workspace_lf: Allocated memory for the operation in vector_implicit_euler() and newton_solver(). Must be of size \f$n\cdot(5+2n)\cdot\text{sizeof}(\text{double})\f$.
 * @param[in] workspace_d: Allocated memory for the DGESV routine. Used for storing permutations. Must be \f$n\cdot\text{sizeof}(\text{int})\f$.
 * @param[in] max_iterations: Maximal number of iterations for each call to the Newton solver.
 * @param[in] tolerance: The desired tolerance for the infinity norm used in the Newton solver.
 * @param[in] f_func: functiontype() pointer to the drift term.
 * @param[in] g_func: functiontype() pointer to diffusion term.
 * @param[in] J_func: functiontype() pointer to the Jacobian of the drift term.
 * @param[in] pu: Pointer to control parameters.
 * @param[in] pd: Pointer to the disturbances.
 * @param[in] pP: Vopid pointer to an arbitrary parameter input.
 * @param[in] px0: Pointer to the initial condition. Must be of size \f$n\cdot\text{sizeof}(\text{double})\f$.
 * 
 * @author Anton Rydahl
 * 
 * @date 20th of November 2020
 * 
 */

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
);

/**
 * The Newton solver implements a root finding procedure. It solves the the matrix system 
 * \f$x_{n+1}-f(x_{n+1})(t_{n+1}-t_{n}) -\psi_n= 0\f$ where \f$\psi_k\f$ consists of the solution in the
 * previous step and the diffusion term in the current step: \f$\psi_n=x_n+g(x_n)d\omega_n\f$.
 * 
 * @param[in] f_func: functiontype() pointer to the drift term.
 * @param[in] J_func: functiontype() pointer to the Jacobian of the drift term.
 * @param[in] max_iterations: Maximal number of iterations used if convergence is not obtained.
 * @param[in] tolerance: The desired tolerance for the infinity norm used in the Newton solver.
 * @param[in] n: The dimension of the system. 
 * @param[in] dt: The size of the temporal step.
 * @param[in] pt: Pointer to the temporal solution. Only used as input parameter to f_func and J_func. Must at least be of size \f$1\cdot\text{sizeof}(\text{double})\f$.
 * @param[in,out] pt: Pointer to the spatial solution. On input this must contain the initial guess for \f$x_{n+1}\f$. 
 * After execution it will contain the final guess for \f$x_{n+1}\f$. Must be of size \f$n\cdot\text{sizeof}(\text{double})\f$.
 * @param[in] ppsi: Should contain \f$\psi_n=x_n+g(x_n)d\omega_n\f$. Must be of size \f$n\cdot\text{sizeof}(\text{double})\f$.
 * @param[in] workspace_lf: Allocated memory of doubles used for storing the residuals. Must be of size \f$n\cdot(2+2\cdot n)\cdot\text{sizeof}(\text{double})\f$.
 * @param[in] workspace_d: Allocated memory of integers used  for storing the row permutation indexes in DGESV. Must be of size \f$n\cdot\text{sizeof}(\text{int})\f$.
 * @param[in] pu: Pointer to control parameters.
 * @param[in] pd: Pointer to the disturbances.
 * @param[in] pP: Vopid pointer to an arbitrary parameter input.
 * 
 * @author Anton Rydahl
 * 
 * @date 20th of November 2020
 * 
 */

void newton_solver(
    functiontype f_func,
    functiontype J_func,
    int max_iterations,
    double tolerance,
    int n, // dimension of x
    double dt,
    double *pt,
    double *px, // input and output
    double *ppsi,
    double *workspace_lf, // double workspace
    int *workspace_d, // integer workspace
    double *pu,
    double *pd,
    void *pP
);

/**
 * Implementation of the implicit-explicit Euler method. This numerical scheme approximates 
 * \f$dx(t) = f\big(x(t)\big)dt+g\big(x(t)\big)d\omega(t)\f$ by 
 * \f$x_{n+1} =x_{n} f(t_{n},x_{n})(t_{n+1}-t_{n})+g(t_{n},x_{n})d\omega_{n}\f$.
 * This implementation differs from vector_implicit_euler() in the sense that it only stores the final step in the solution.
 * This ix way more space efficient
 * 
 * @param[in] N: The number of time steps (excluding the initial condition).
 * @param[in] n: The dimension of \f$x(t)\f$.
 * @param[in] NS: The number of simulations.
 * @param[in] pt: Pointer to the temporal solution.  Must be of size \f$(N+1)\cdot\text{sizeof}(\text{double})\f$.
 * @param[out] px: After operation this array will contain the spatial solution. Must be size \f$n\cdot (N+1)\cdot NS\cdot\text{sizeof}(\text{double})\f$.
 * @param[in] pdW: White noise. Must be size \f$n\cdot (N+1)\cdot NS \cdot\text{sizeof}(\text{double})\f$.
 * @param[in] workspace_lf: Allocated memory for the operation in vector_implicit_euler() and newton_solver(). Must be of size \f$n\cdot(7+2n)\cdot\text{sizeof}(\text{double})\f$.
 * @param[in] workspace_d: Allocated memory for the DGESV routine. Used for storing permutations. Must be \f$n\cdot\text{sizeof}(\text{int})\f$.
 * @param[in] max_iterations: Maximal number of iterations for each call to the Newton solver.
 * @param[in] tolerance: The desired tolerance for the infinity norm used in the Newton solver.
 * @param[in] f_func: functiontype() pointer to the drift term.
 * @param[in] g_func: functiontype() pointer to diffusion term.
 * @param[in] J_func: functiontype() pointer to the Jacobian of the drift term.
 * @param[in] pu: Pointer to control parameters.
 * @param[in] pd: Pointer to the disturbances.
 * @param[in] pP: Void pointer to an arbitrary parameter input.
 * @param[in] px0: Pointer to the initial condition. Must be of size \f$n\cdot\text{sizeof}(\text{double})\f$.
 * 
 * @author Anton Rydahl
 * 
 * @date 10th of December 2020
 * 
 */


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
);

#endif
