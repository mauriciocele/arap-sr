#if !defined(__SVD__H__)
#define __SVD__H__

//Solves A·X = B for a vector X, where A is specified by the arrays u[1..m][1..n], w[1..n],
//v[1..n][1..n] as returned by svdcmp. m and n are the dimensions of a, and will be equal for
//square matrices. b[1..m] is the input right-hand side. x[1..n] is the output solution vector.
//No input quantities are destroyed, so the routine may be called sequentially with different b’s.
void svbksb(double **u, double *w, double **v, int m, int n, double *b, double *x);

//Given a matrix A[1..m][1..n], this routine computes its singular value decomposition, A =
//U·W·V^T. Thematrix U replaces A on output. The diagonal matrix of singular values W is output
//as a vector w[1..n]. Thematrix V (not the transpose V^T ) is output as v[1..n][1..n].
void svdcmp(double **a, int m, int n, double *w, double **v);

void svdcmp3x3(double a[4][4], int m, int n, double w[], double v[4][4]);

#endif
