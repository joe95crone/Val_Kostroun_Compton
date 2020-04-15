#include <math.h>
#include <stdlib.h>
//#include <iomanip.h>
//#include <stdio.h>
#include "..\Header\matrix2.h"



/******************************************************************/
void init_matrix(double **A, int nrow, int ncol)
/*
Takes a matrix A[0..nrow-1][0...ncol-1] and initialized all elements to zero
*/
{
int i, j;

for (i=0;i<nrow;i++)
{
    for (j=0;j<ncol;j++)
    {
        A[i][j] = 0.0;
    }
}

}//init_matrix
/******************************************************************/
double **dtranspose(double **A, int nrow, int ncol)
/* Given matrix A[0..nrow][0..ncol], returns the transpose of A */
{
int i,j;
double **B;
B = dmatrix(nrow,ncol);
for (i = 0;i<nrow;i++)
{
    for (j = 0;j<ncol;j++)
    {
        B[i][j] = A[j][i];
    }
}

return(B);
}//dtranspose

/******************************************************************/

double **dmatrix_mult(double ** A, double ** B, int n_row_a, int n_col_a_n_row_b, int n_col_b)
/* Given matrix A(n_row_a, n_col_a) and matrix B(n_row_b, n_col_b) where n_col_a = n_row_b = n_col_a_n_row_b
this function returns a new matrix or vector C(n_row_a, n_col_b) = A*B */
{
int i, j, k;
double **C;
C = dmatrix(n_row_a,n_col_b);
init_matrix(C,n_row_a,n_col_b);
for (i=0; i<n_row_a; i++)
{
    for (j = 0; j<n_col_b;j++)
    {
        for (k=0;k<n_col_a_n_row_b;k++)
        {
            C[i][j] = C[i][j] + A[i][k]*B[k][j];
        }
    }
}
return(C);
}//dmatrix_mult

/*******************************************************************************************************/
double** copy_matrix(double**A, int nrow, int ncol)
{
int i,j;
double **B;
B = dmatrix(nrow,ncol);
for (i = 0;i<nrow;i++)
{
    for (j = 0;j<ncol;j++)
    {
        B[i][j] = A[i][j];
    }
}

return(B);

}


/******************************************************************/


void choldc(double** a, int n)
/*
Given a positive-definite symmetric matrix a[0..n-1][0..n-1], this routine constructs the Cholesky
decomposition, A = L*L^T. Taken from numerical recipies, slightly modified.
*/
{
    //void nrerror(char error_text[]);
    int i,j,k;
    double sum;
    double *p;
    p= new double[n];
    //p = malloc(n*sizeof(double));

    for (i=0;i<n;i++) {
        for (j=i;j<n;j++) {
            for (sum=a[i][j],k=i-1;k>=0;k--) sum -= a[i][k]*a[j][k];
            if (i == j) {
            //  if (sum <= 0.0)
            //      nrerror("choldc failed");
                p[i]=sqrt(sum);
            } else a[j][i]=sum/p[i];
        }
    }
    for (i=0;i<n;i++)
    {
        a[i][i]= p[i];
        for (j=i+1;j<n;j++)
        {
            a[i][j] =0;
        }
    }
    free(p);
}


/*********************************************************************************/

double **dmatrix(long dim_r, long dim_c)
/* allocate a double matrix with subscript range m[0..dim_r-1][0..dim_c-1] */
{
    long nrl = 0;
    long ncl = 0;
    long nrh = dim_r -1;
    long nch = dim_c - 1;
    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    double **m;

    /* allocate pointers to rows */
    m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
    //if (!m) nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
    //if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

    /* return pointer to array of pointers to rows */

    return m;
}
/***********************************************************************************/
void free_dmatrix(double **m, long dim_r, long dim_c)
/* free a double matrix allocated by dmatrix() */
{
	long nrl = 0;
    long ncl = 0;
    long nrh = dim_r -1;
    long nch = dim_c - 1;
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}


/***********************************************************************************/

/* void display_matrix(double **A, int nrow, int ncol)
{
    int ii, jj;
    for(ii=0;ii<nrow;ii++)
    {
        for(jj=0;jj<ncol;jj++)
        {
            cout<<setw(10)<<A[ii][jj];
        }
        cout<<endl;
    }
    cout<<endl;
}*/

/*************************************************************************************/
int mat_ind (int i, int j, int n_rows)
/* Takes in row index, i, and col index, j, of a 2D matrix and returns the corresponding
array index assuming n_rows rows.
*/
{
	int temp;

	temp = i + (j*n_rows);
	return(temp);

}
