#ifndef RAND_H
#define RAND_H 1




/**************Definitions for ran2 function*************************/
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
/**********************************************************************/




double ran2(long *);
/*
Long period (>2x10^18) random number generator of L'Ecuyer with Bays-Durham shuffle
and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of 
the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter idum
between successive deviates in a sequence.  RNMX should approximate the largest floating value that is less
than 1.
*/


double randn2(long *);
//returns a random number taken from a normal distribution with zero mean and unity standard deviation.
// If init ~= 0, then the random number generator will be initizlized with a call to time. Otherwise, it is not
// initiazlized.


long Get_ran2_init();
/*
Returns a negative long integer obtained from the time function. The address of
(pointer to) the returned value should be used to initialize calls to ran2 or
randn2.
*/


int randnCovA(double **, double[], int, long*);
/*
Takes square, positive definite covariance matarix CovA, for a Gaussian distributed ndim dimensional state
vector, and caclulates a single instance of a random draw of the state vector from the distriubtion described
by CovA.
*/




#endif
