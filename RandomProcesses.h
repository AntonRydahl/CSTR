/// @file RandomProcesses.h

#ifndef SDE_LIBRARY
#define SDE_LIBRARY

#include <math.h>

///@cond IGNORE_IN_DOCUMENTATION

extern void dtrmm_(char * side, char * uplo, char * transa, char * diag,
              int * m, int * n, double * alpha, double * A, int * lda,
              double * B, int * ldb);

extern void dpotrf_(char * uplo,int * N, double * A,int * lda, int * info);

///@endcond

/**
 * This method generates stochastic variables following a distribution
 * where the mean vector is the zero vector and the covariance matrix 
 * can be describes as dt * I. Here I is the identity matrix and dt is
 * the temporal derivative
 * Hence the cholesky factorization of this matrix is sqrt(dt) * I 
 * and one does not have to provide the covariance matrix as input.
 * This method uses column major storage.
 * 
 * @param[out] pdW: White noise. Must be of size \f$n_\omega\cdot N\cdot NS\cdot\text{sizeof}(\text{double})\f$.
 * @param[in] pworkspace: Workspace for Mercenne Twister. Must be of size \f$624\cdot\text{sizeof}(\text{unsigned long})\f$.
 * @param[in] FT: Final time.
 * @param[in] N: The number of time intervals. The number of columns in pT, pdW and W. 
 * @param[in] nw: The dimension of the Wiener process. The number of rows in pdW and W.
 * @param[in] NS: The number of simulations. The number of two dimensional matrices in pdW and W.
 * 
 * @author Anton Rydahl
 * @date 18th of September 2020
 * 
 */

void scalar_wiener_process(
    double *pdW,
    unsigned long *pworkspace,
    double FT,
    int N,
    int nw,
    int NS
);

/**
 * This method accumulates the random noise generated with the above methods into a standard Wiener process.
 * The standard Wiener process is the cummulative sum of the white noise over time, and hence this function computes
 * the row-wise cummulative sum of the matrix pdW and stores the result in W. Note that the first column will contain only 
 * zeros.
 * 
 * @param[in,out] pW: Standard Wiener process. Must be of size \f$n_\omega\cdot(N+1)\cdot NS\cdot\text{sizeof}(\text{double})\f$.
 * @param[in] pdW: White noise. Must be of size \f$n_\omega\cdot N\cdot NS\cdot\text{sizeof}(\text{double})\f$.
 * @param[in] N: The number of time intervals. The number of columns in pT, pdW and W. 
 * @param[in] nw: The dimension of the Wiener process. The number of rows in pdW and W.
 * @param[in] NS: The number of simulations. The number of two dimensional matrices in pdW and W.
 * 
 * @author Created by Anton Rydahl
 * @date 27th of Sepetember 2020
 */

void cumsum(
    double *pW,
    double *pdW,
    int N,
    int nw,
    int NS
);

/**
 * This method generated an equidistant vector of length \f$N+1\f$ consisting of \f$\{t_0,t_0+dt,\ldots,t_0+T\}\f$ where \f$dt=\frac{T}{N}\f$.
 * 
 * @param[in,out] pT: Temporal solution. Must be of size \f$(N+1)\cdot\text{sizeof}(\text{double})\f$.
 * @param[in] t0: Initial time.
 * @param[in] FT: Final time.
 * @param[in] N: The number of time intervals. The number of columns in pT, pdW and W. 

 * @author Created by Anton Rydahl
 * @date 1st of November 2020
 */

void linspace(
    double *pT,
    double t0,
    double FT,
    int N
);

#endif
