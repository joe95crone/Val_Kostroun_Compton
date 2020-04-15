#include <math.h>
//#include <fstream.h>
//#include <iomanip.h>
//#include <string.h>
#include <stdio.h>

////////////////////////
#include <stdlib.h>
#include <time.h>

#include "..\Header\random2.h"
#include "..\Header\matrix2.h"
///////////////////////





/*******************************************************************************************************/

int randnCovA(double **CovA, double Y[], int ndim, long* init)
/*
Takes square, positive definite covariance matarix CovA, for a Gaussian distributed ndim dimensional state
vector, and caclulates a single instance of a random draw of the state vector from the distriubtion described
by CovA.
*/
{
double ** B;
int i,j;

double *X;
X= new double[ndim];


//Initialize Y vector
for (i = 0; i < ndim; i++)
{
    Y[i] = 0.0;
}

B = copy_matrix(CovA,ndim,ndim);
//Peform Cholesky decomposition to get factor A=B*B^T
choldc(B,ndim);

//Get random draw on state vector with covariance matrix = identity matrix
for (i=0;i<ndim;i++)
{
    X[i]=randn2(init);
}

//Multiply by B to get properly transformed random variable.
for (i=0;i<ndim;i++)
{
    for (j=0;j<ndim;j++)
    {
        Y[i] = Y[i] + B[i][j]*X[j];
    }//for j

}//for i


free (X);

return(0);

}//randnCovA

/*******************************************************************************************************/

double randn2(long* init)
//returns a random number taken from a normal distribution with zero mean and unity standard deviation.
// If init ~= 0, then the random number generator will be initizlized with a call to time. Otherwise, it is not
// initiazlized.
{

    double num1, num2, ans1;
    //double ans2;
    double pi2 = 4.*acos(0.0);




/*  if(num1 < ZERO)
    {
        num1 = rand1(0)*(1./double(RAND_MAX));
        cout<<num1<<endl;
    }*/


    num1 = ran2(init);
    num2 = ran2(init);

    ans1 = sqrt(-2.*log(num1))*cos(pi2*num2);


    //ans2 = sqrt(-2.*log(num1))*sin(pi2*num2);

    return(ans1);


}//rand1

/*******************************************************************************************************/



double ran2(long *idum)
{
/*
Long period (>2x10^18) random number generator of L'Ecuyer with Bays-Durham shuffle
and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter idum
between successive deviates in a sequence.  RNMX should approximate the largest floating value that is less
than 1.
*/

    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    double temp;

    if (*idum <= 0) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        idum2=(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}


/******************************************************************/

long Get_ran2_init()
/*
Returns a negative long integer obtained from the time function. The address of
(pointer to) the returned value should be used to initialize calls to ran2 or
randn2.
*/
{
    return(-time(0));
}
