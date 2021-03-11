/// @file MersenneTwister.h

#ifndef DOUBLE_PRECISION_RANDOM_NUMBERS
#define DOUBLE_PRECISION_RANDOM_NUMBERS

/**
 * This function will set the inital value (seed) used in both mersenne_twister() and in the
 * linear_congruential(). Note that these functions continuously update the seed such that they
 * will generate a new distribution of random numbers every time they are invoked. Hence the seed must be set before every call
 * to the uniform random number generators, if one wants to control the seed.
 * 
 * @param[in] x: This will be set as the new seed.
 * 
 * @author Anton Rydahl
 * @date 30th of August 2020
*/

void set_seed(
    unsigned int x
);

/**
 * This is double precision implementation of the Linear Congruential Generator which can be used for both singe, double and extended precision.
 * This pseudo random number generator generates a sample taken from the uniform distribution on the unit interval \f$(0,1)\f$. 
 * The statistical properties of the sample generated are in general not as good as for Mersenne Twister , but it is a way faster method for
 * generating pseudo random numbers. It has been chosen to implement the method with modulus \f$m=2^{32}\f$ together with the increment 
 * \f$c=2531011\f$ and the multiplier \f$a=1103515245\f$. These are the same as in the Microsoft Visual C++ standard from 1993 the LCG and guarantees
 * that the statistical properties of the distribution does not depend on the seed.
 * 
 * @param[out] parr: Pointer to memory of length \f$N\cdot \text{sizeof}(double)\f$.
 * @param[in] N: Number of random variables to be generated.
 * 
 * @author Anton Rydahl
 * @date 30th of August 2020
*/

void linear_congruential(
    double *parr,
    int N
);

/**
 * This is the 32 bit implementation of the Mersenne Twister algorithm similar to the
 * <a href="http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/980409/mt19937int.out">original 1998 edition</a>. 
 * It can be used to generate a sample of pseudo random numbers drawn from the uniform distribution on the the open unit domain \f$(0,1)\f$.
 * There are no big changes in the algorithm compared to the 1998 edition (except for the output domain), but the implementation has been changed to avoid using static arrays
 *  and to provide a more userfriendly interface that can be used for generating floating points of single, double and extended precision. 
 * 
 * The array "pworkspace" is used to store the generator of the algorithm and can be reused in further calls to this routine. 
 * 
 * @param[out] parr: Pointer to memory of length \f$N\cdot \text{sizeof}(double)\f$.
 * @param[in] pgenerator: Pointer to memory of length \f$624\cdot\text{sizeof}(\text{unsigned long})\f$. This is used for storing the generator.
 * @param[in] N: Number of random variables to be generated.
 * 
 * @author Anton Rydahl
 * @date 30th of August 2020
*/

void mersenne_twister(
    double *parr, 
    unsigned long *pgenerator, 
    int N
);

/**
 * This is the implementation of the Box-Muller transformation which transforms uniformly distributed random variables into normally
 * distributed random variables in sets of two. Note that if an uneven number of variables is inputted (if 2 divides N-1), the function 
 * will reuse the first uniform random number in order to trasform the last variable.
 * 
 * Also, please note that "parr" has to be initialized with either linear_congruential() or mersenne_twister() before 
 * being fed to this method. It simply assumes this input already holds uniformly distributed random variables.
 * 
 * @param[in,out] parr: Pointer to memory of length \f$N\cdot \text{sizeof}(double)\f$.
 * @param[in] N: Number of random variables held by parr.
 * @param[in] mu: The desired scalar mean of the distribution.
 * @param[in] sigma: The desired scalar standard deviation of the distribution.
 * 
 * @author Anton Rydahl
 * @date 2nd of September 2020
*/

void box_muller(
    double *parr,
    int N,
    long double mu, 
    long double sigma
);

/**
 * A double precision implementation of the rejection based ratio of unitforms normal transformation. 
 * Since this is a rejection based method, it does not always need exactly the same amount of uniform random numbers. 
 * Therefore you need to provide an additional workspace for the algorithm. "pworksapce" is used for storing temporary
 * uniform random numbers which are then transformed and stored in "parr". 
 * 
 * @param[out] parr: Pointer to memory of length \f$N\cdot \text{sizeof}(double)\f$.
 * @param[in] N: Number of random variables held by parr.
 * @param[in] mu: The desired scalar mean of the distribution.
 * @param[in] sigma: The desired scalar standard deviation of the distribution.
 * @param[in] pgenerator: Pointer to memory of length \f$624\cdot\text{sizeof}(\text{unsigned long})\f$. For use in Mersenne Twister.
 * @param[in] pworkspace: Workspace of size \f$N\cdot \text{sizeof}(double)\f$.
 * 
 * @author Anton Rydahl
 * @date 2nd of September 2020
*/

void ratio_of_uniforms(
    double *parr,
    int N,
    long double mu, 
    long double sigma,
    unsigned long *pgenerator,
    double *pworkspace
);

/**
 * A double precision implementation of the rejection based normal transformation called normal polar transformation. 
 * Since this is a rejection based method, it cannot be predicted how many uniform random numbers the method will need.
 * Therefore the method generates uniform random numbers with mersenne_twister(). This means that you need to provide 
 * an additional workspace for the algorithm. "pworksapce" is used for storing temporary
 * uniform random numbers which are then transformed and stored in "parr". 
 * 
 * @param[out] parr: Pointer to memory of length \f$N\cdot \text{sizeof}(double)\f$.
 * @param[in] N: Number of random variables held by parr.
 * @param[in] mu: The desired scalar mean of the distribution.
 * @param[in] sigma: The desired scalar standard deviation of the distribution.
 * @param[in] pgenerator: Pointer to memory of length \f$624\cdot\text{sizeof}(\text{unsigned long})\f$. For use in Mersenne Twister.
 * @param[in] pworkspace: Workspace of size \f$N\cdot \text{sizeof}(double)\f$.
 * 
 * @author Anton Rydahl
 * @date 2nd of September 2020
*/

void normal_polar(
    double * parr,
    int N,
    long double mu, 
    long double sigma,
    unsigned long *pgenerator,
    double *pworkspace
);

/**
 * This routine generates double precision uniformly distributed random variables on the unit interval. 
 * It invokes the Mersenne Twister implementation to generate a uniform sample.
 * 
 * @param[out] parr: Pointer to memory of length \f$N\cdot\text{sizeof}(\text{double})\f$.
 * @param[in] pworkspace: Pointer to memory of length \f$624\cdot\text{sizeof}(\text{unsigned long})\f$.
 * @param[in] N: Number of random variables to be generated.
 * 
 * @author Anton Rydahl
 * @date 9th of Sepember 2020
*/

void d_rand_uniform(
    double *parr, 
    unsigned long *pworkspace, 
    int N
);

/**
 * Double precision standard normally distributed random variables. The routine uses the
 * Mersenne Twister algorithm to generate a uniform distribution of random variables
 * which is transformed into a normal distribution using Box Mueller transformation.
 * 
 * @param[out] parr: Pointer to memory of length \f$N\cdot\text{sizeof}(\text{double})\f$.
 * @param[in] pworkspace: Pointer to memory of length \f$624\cdot\text{sizeof}(\text{unsigned long})\f$.
 * @param[in] N: Number of random variables to be generated.
 * 
 * @author Anton Rydahl 
 * @date 9th of Sepember 2020
*/

void d_rand_standard_normal(
    double *parr, 
    unsigned long *pworkspace, 
    int N
);

/**
 * Double precision normally distributed random variables. This routine uses the
 * Mersenne Twister algorithm to generate a uniform distribution of random variables
 * which is transformed into a normal distribution using Box Mueller transformation.
 * 
 * @param[out] parr: Pointer to memory of length \f$N\cdot\text{sizeof}(\text{double})\f$.
 * @param[in] pworkspace: Pointer to memory of length \f$624\cdot\text{sizeof}(\text{unsigned long})\f$.
 * @param[in] N: Number of random variables to be generated.
 * @param[in] mu: Mean of distribution.
 * @param[in] sigma: Standard deviation of distribution.
 * 
 * @author Anton Rydahl
 * @date 9th of Sepember 2020
*/

void d_rand_normal(
    double *parr, 
    unsigned long *pworkspace, 
    int N,
    long double mu, 
    long double sigma
);

#endif
