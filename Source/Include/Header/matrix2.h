#ifndef MATR_H
#define MATR_H 1

#define NR_END 1
#define FREE_ARG char*



double **dmatrix(long, long);
/* allocate a double matrix with subscript range m[0..dim_r-1][0..dim_c-1] */

void free_dmatrix(double **m, long, long);

void init_matrix(double **, int, int);
/* 
Takes a matrix A[0..nrow-1][0...ncol-1] and initialized all elements to zero
*/


double **dtranspose(double **, int, int);
/* Given matrix A[0..nrow][0..ncol], returns the transpose of A */

double** copy_matrix(double**, int, int);


double **dmatrix_mult(double **, double **, int, int, int);
/* Given matrix A(n_row_a, n_col_a) and matrix B(n_row_b, n_col_b) where n_col_a = n_row_b = n_col_a_n_row_b
this function returns a new matrix or vector C(n_row_a, n_col_b) = A*B */

void choldc(double **, int);
/*
Given a positive-definite symmetric matrix a[0..n-1][0..n-1], this routine constructs the Cholesky
decomposition, A = L*L^T. Taken from numerical recipies, slightly modified.
*/

//void display_matrix(double **, int, int);


int mat_ind (int, int, int);
/* Takes in row index, i, and col index, j, of a 2D matrix and returns the corresponding
array index assuming n_rows rows.
*/



#endif 