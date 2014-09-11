/*--------------------------------------------------------------------------*/
/** \file mat_sym2.h   

   \brief Basic operations with symmetric 2 * 2 matrices. 

   The matrices are stored as a field of 6 doubles such that a field 
   \f$ \tilde{A}=(a[0], a[1], a[2]) \f$ stands for the matrix
   \f$ A = 
     \left(
       \begin{array}{cc}
         a[0] & a[2] \\
	 a[2] & a[1] \\
       \end{array}
     \right)
   \f$.
   In the following documentation, we refer to \f$\tilde{A}\f$ as 
   vector format and to \f$A\f$ as matrix format. 

   Maintainer: Stephan Didas <didas@mia.uni-saarland.de>                    */
/*--------------------------------------------------------------------------*/

#ifndef __MAT_SYM2_H__
#define __MAT_SYM2_H__

#include <stdlib.h>

#include "base/memory.h"
#include "matrix/matrix.h"

/*--------------------------------------------------------------------------*/

/** \brief Adds two symmetric 2*2 matrices and stores the sum in the first
           one. 
 */
void mat_sym2_add(double* a,  /**< the first summand  */
		  double* b)  /**< the second summand */
{
    a[0] = a[0] + b[0];
    a[1] = a[1] + b[1];
    a[2] = a[2] + b[2];

    return;
} /* mat_sym2_add */

/*--------------------------------------------------------------------------*/

/** \brief Subtracts one symmetric 2*2 matrix from another and stores the 
           result in the first one.
 */
void mat_sym2_sub(double* a,  /**< minuend and result (in- and output) */
		  double* b)  /**< subtrahend                          */
{
    a[0] = a[0] - b[0];
    a[1] = a[1] - b[1];
    a[2] = a[2] - b[2];

    return;
} /* mat_sym2_sub */

/*--------------------------------------------------------------------------*/

/** \brief Adds two symmetric 2*2 matrices and stores the sum in a third
           one. 
 */
void mat_sym2_add_get(double* a,   /**< the first summand, input   */
		      double* b,   /**< the second summand, input  */
		      double* res) /**< the sum of a and b, output */
{
    res[0] = a[0] + b[0];
    res[1] = a[1] + b[1];
    res[2] = a[2] + b[2];

    return;
} /* mat_sym2_add_get */

/*--------------------------------------------------------------------------*/

/** \brief Subtracts one symmetric 2*2 matrix from another and stores the 
           result in a third one.
 */
void mat_sym2_sub_get(double* a,   /**< minuend and result, input */
		      double* b,   /**< subtrahend, input         */
		      double *res) /**< the difference, output    */
{
    res[0] = a[0] - b[0];
    res[1] = a[1] - b[1];
    res[2] = a[2] - b[2];

    return;
} /* mat_sym2_sub_get */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the mean value of two symmetric 2*2 matrices.
 */
void mat_sym2_mean_get(double* a,   /**< minuend and result, input */
		       double* b,   /**< subtrahend, input         */
		       double *res) /**< the difference, output    */
{
    double one_half = 0.5;

    res[0] = ( a[0] + b[0] ) * one_half;
    res[1] = ( a[1] + b[1] ) * one_half;
    res[2] = ( a[2] + b[2] ) * one_half;

    return;
} /* mat_sym2_mean_get */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the square of a symmetric 2*2 matrix.
 */
void mat_sym2_sqr(double* arg,  /**< the argument */
		  double* res)  /**< the result   */
{
    res[0] = arg[0]*arg[0] + arg[2]*arg[2];
    res[1] = arg[1]*arg[1] + arg[2]*arg[2];
    res[2] = arg[0]*arg[2] + arg[1]*arg[2];

    return;
} /* mat_sym2_sqr */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the symmetric product of two symmetric 2*2 matrices.

    The symmetric product of \f$A\f$ and \f$B\f$ is defined here as
    \f$A B A\f$.

    Note: In most of the cases, it will be necessary to calculate the 
          square root of A first.
 */
void mat_sym2_symm_prod(double* a,   /**< matrix A                    */
			double* b,   /**< matrix B                    */
			double* res) /**< the symmetrix product A*B*A */
{
    res[0] = ( a[0]*b[0] + a[2]*b[2] ) * a[0]
	+ ( a[0]*b[2] + a[2]*b[1] ) * a[2];
    res[1] = ( a[2]*b[0] + a[1]*b[2] ) * a[2]
	+ ( a[2]*b[2] + a[1]*b[1] ) * a[1];
    res[2] = ( a[2]*b[0] + a[1]*b[2] ) * a[0]
	+ ( a[2]*b[2] + a[1]*b[1] ) * a[2];

    return;
} /* mat_sym2_symm_prod */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the Jordan product of two symmetric 2*2 
           matrices.

    This product is defined as \f$\frac{1}{2} (AB + BA)\f$.
 */
void mat_sym2_jordan_prod(double *a,   /**< matrix A */
			  double *b,   /**< matrix B */
			  double *res) /**< the product */
{
    res[0] = a[0]*b[0] + a[2]*b[2];
    res[1] = a[1]*b[1] + a[2]*b[2];
    res[2] = 0.5 * ( a[0]*b[2] + a[2]*b[1] + a[2]*b[0] + a[1]*b[2] );
    
    return;
} /* mat_sym2_jordan_prod */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the Schur product of two symmetric 2*2 matrices.

    The Schur product is the componentwise product of the two matrices. 
 */
void mat_sym2_schur_prod(double *a,   /**< matrix A */
			 double *b,   /**< matrix B */
			 double *res) /**< the product */
{
    res[0] = a[0]*b[0];
    res[1] = a[1]*b[1];
    res[2] = a[2]*b[2];
    
    return;
} /* mat_sym2_schur_prod */

/*--------------------------------------------------------------------------*/

/** \brief Prints a symmetric 2*2 matrix out on screen. 
 */
void mat_sym2_print(double* a)  /**< the matrix */
{
    printf("  ( %18.10e  %18.10e )\n", a[0], a[2] );
    printf("  ( %18.10e  %18.10e )\n", a[2], a[1] );

    return;
} /* mat_sym2_print */

/*--------------------------------------------------------------------------*/

/** \brief Copies a symmetric 2*2 matrix into another one. 
 */
void mat_sym2_copy(double *source,  /**< source matrix      */
		   double *dest)    /**< destination matrix */
{
    dest[0] = source[0];
    dest[1] = source[1];
    dest[2] = source[2];

    return;
} /* mat_sym2_copy */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the Frobenius norm of a symmetric 2*2 matrix.
 */
double mat_sym2_frobnorm(double *a)  /**< the matrix */
{
    return sqrt( a[0]*a[0] + a[1]*a[1] + 2.0*a[2]*a[2] );
} /* mat_sym2_frobnorm */

/*--------------------------------------------------------------------------*/

/** \brief Calculates a + lambda b for two symmetric 2*2 matrices and a 
           scalar lambda and stores the result in a third matrix r. 
 */
void mat_sym2_a_add_lb_get(double* a,      /**< the first summand  */
			   double  lambda, /**< the scalar factor  */ 
			   double* b,      /**< the second summand */
			   double* r)      /**< the result         */
{
    r[0] = a[0] + lambda * b[0];
    r[1] = a[1] + lambda * b[1];
    r[2] = a[2] + lambda * b[2];

    return;
} /* mat_sym2_a_add_lb_get */

/*--------------------------------------------------------------------------*/
/** \brief Calculates a + lambda b for two symmetric 2*2 matrices and a 
           scalar lambda and stores the result in the first matrix a. 
 */
void mat_sym2_a_add_lb(double* a,      /**< the first summand  */
		       double  lambda, /**< the scalar factor  */ 
		       double* b)      /**< the second summand */
{
    a[0] = a[0] + lambda * b[0];
    a[1] = a[1] + lambda * b[1];
    a[2] = a[2] + lambda * b[2];

    return;
}

/*--------------------------------------------------------------------------*/

/** \brief Converts a symmetric 2*2 matrix from vector format to matrix 
           format.

    This method is necessary to use the PA_trans method from matrix.c.
    For the explanation of vector and matrix format see also the head of 
    documentation of the file mat_sym3.h.
 */
void mat_sym2_to_matrix(double *in,  /**< the matrix in vector format   
                                          unchanged (input)             */
			double **b)  /**< the matrix in matrix format, 
					  output                        */
{
    /* copy input to b */
    b[1][1] = in[0];
    b[1][2] = in[2];
    b[2][2] = in[1];

    b[2][1] = in[2];

    return;
} /* mat_sym2_to_matrix */

/*--------------------------------------------------------------------------*/

/** \brief Converts a symmetric 2*2 matrix from matrix format to vector
           format. 

    This method is necessary to use the PA_backtrans function from matrix.c.
    For the explanation of vector and matrix format see also the head of 
    documentation of the file mat_sym2.h.
 */
void mat_sym2_from_matrix(double **b,   /**< the matrix in matrix format */
			  double *out)  /**< the matrix in vector format */
{
    /* copy result to out */
    out[0] = b[1][1];
    out[1] = b[2][2];
    out[2] = b[1][2];

    return;
} /* mat_sym2_from_matrix */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the square of the negative part of a symmetric 2*2
           matrix. 

    Computes a principal axis transformation and sets all nonnegative 
    eigenvalues to zero. The negative ones are squared, and the resulting
    matrix is calculated from the new eigenvalues and the eigenvectors of 
    the input matrix. 
 */
void mat_sym2_min_zero_sq(double* in,   /**< input matrix, unchanged */
			  double* out)  /**< result, changed         */
{
    double** b;      /* storage for matrix (for PA method) */
    double** v;      /* the eigenvectors */
    double*  lambda; /* the eigenvalues  */

    double delta = 10e-14;
    double eps   = 10e-08;

    matrix_alloc(&b, 4, 4);
    matrix_alloc(&v, 4, 4);
    vector_alloc(&lambda, 4);

    mat_sym2_to_matrix(in, b);

    PA_trans(b, 2, delta, eps, lambda, v);

    if(lambda[1] < 0) {
	lambda[1] = lambda[1] * lambda[1];
    } else {
	lambda[1] = 0.0;
    }

    if(lambda[2] < 0) {
	lambda[2] = lambda[2] * lambda[2];
    } else {
	lambda[2] = 0.0;
    }

    PA_backtrans(lambda, v, 2, b);

    mat_sym2_from_matrix(b, out);

    matrix_disalloc(b);
    matrix_disalloc(v);
    vector_disalloc(lambda);

    return;
} /* mat_sym2_min_zero_sq */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the square root of a symmetric 2*2 matrix. 

    Computes a principal axis transformation and takes the square root of 
    all eigenvalues. Then resulting matrix is calculated from the new 
    eigenvalues and the eigenvectors of the input matrix. 
     
 */
void mat_sym2_sqrt(double* in,  /**< input matrix, unchanged */
		   double* out) /**< result, changed         */
{
    double** b;      /* storage for matrix (for PA method) */
    double** v;      /* the eigenvectors */
    double*  lambda; /* the eigenvalues  */

    double delta = 10e-10;
    double eps   = 10e-06;

    matrix_alloc(&b, 4, 4);
    matrix_alloc(&v, 4, 4);
    vector_alloc(&lambda, 4);

    mat_sym2_to_matrix(in, b);

    PA_trans(b, 2, delta, eps, lambda, v);

    lambda[1] = sqrt( fabs(lambda[1]) );
    lambda[2] = sqrt( fabs(lambda[2]) );

    PA_backtrans(lambda, v, 2, b);

    mat_sym2_from_matrix(b, out);

    matrix_disalloc(b);
    matrix_disalloc(v);
    vector_disalloc(lambda);

    return;
} /* mat_sym2_sqrt */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the reciprocal of the square root of a symmetric 
           2*2 matrix. 

    Computes a principal axis transformation and takes the square root of 
    all eigenvalues. Then resulting matrix is calculated from the new 
    eigenvalues and the eigenvectors of the input matrix. 
 */
void mat_sym2_rec_sqrt(double* in,  /**< input matrix, unchanged */
		       double* out) /**< result, changed         */
{
    double** b;      /* storage for matrix (for PA method) */
    double** v;      /* the eigenvectors */
    double*  lambda; /* the eigenvalues  */

    double delta = 10e-10;
    double eps   = 10e-06;

    matrix_alloc(&b, 4, 4);
    matrix_alloc(&v, 4, 4);
    vector_alloc(&lambda, 4);

    mat_sym2_to_matrix(in, b);

    PA_trans(b, 2, delta, eps, lambda, v);

    if(lambda[1] != 0)
    {
	lambda[1] = 1.0 / sqrt( fabs(lambda[1]) );
    }
    if(lambda[2] != 0)
    {
	lambda[2] = 1.0 / sqrt( fabs(lambda[2]) );
    }

    PA_backtrans(lambda, v, 2, b);

    mat_sym2_from_matrix(b, out);

    matrix_disalloc(b);
    matrix_disalloc(v);
    vector_disalloc(lambda);

    return;
} /* mat_sym2_rec_sqrt */

/*--------------------------------------------------------------------------*/

/** \brief Calculates (fractional) powers of a symmetric 2*2 matrix. 

    Computes a principal axis transformation and takes the power of 
    all eigenvalues. Then resulting matrix is calculated from the new 
    eigenvalues and the eigenvectors of the input matrix. 
 */
void mat_sym2_pow(double* in,  /**< input matrix, unchanged */
		  double  exp, /**< the exponent            */
		  double* out) /**< result, changed         */
{
    double** b;      /* storage for matrix (for PA method) */
    double** v;      /* the eigenvectors */
    double*  lambda; /* the eigenvalues  */

    double delta = 10e-10;
    double eps   = 10e-06;

    matrix_alloc(&b, 4, 4);
    matrix_alloc(&v, 4, 4);
    vector_alloc(&lambda, 4);

    mat_sym2_to_matrix(in, b);

    PA_trans(b, 2, delta, eps, lambda, v);

    lambda[1] = pow( lambda[1], exp );
    lambda[2] = pow( lambda[2], exp );

    PA_backtrans(lambda, v, 2, b);

    mat_sym2_from_matrix(b, out);

    matrix_disalloc(b);
    matrix_disalloc(v);
    vector_disalloc(lambda);

    return;
} /* mat_sym2_pow */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the square of the positive part of a symmetric 2*2 
           matrix. 

    Computes a principal axis transformation and sets all nonpositive     
    eigenvalues to zero. The positive ones are squared, and the resulting 
    matrix is calculated from the new eigenvalues and the eigenvectors of 
    the input matrix. 
 */
void mat_sym2_max_zero_sq(double* in,  /**< input matrix, unchanged */
			  double* out) /**< result, changed         */
{
    double** b;      /* storage for matrix (for PA method) */
    double** v;      /* the eigenvectors */
    double*  lambda; /* the eigenvalues  */

    double delta = 10e-10;
    double eps   = 10e-06;

    matrix_alloc(&b, 4, 4);
    matrix_alloc(&v, 4, 4);
    vector_alloc(&lambda, 4);

    mat_sym2_to_matrix(in, b);

    PA_trans(b, 2, delta, eps, lambda, v);

    if(lambda[1] > 0) {
	lambda[1] = lambda[1] * lambda[1];
    } else {
	lambda[1] = 0.0;
    }

    if(lambda[2] > 0) {
	lambda[2] = lambda[2] * lambda[2];
    } else {
	lambda[2] = 0.0;
    }

    PA_backtrans(lambda, v, 2, b);

    mat_sym2_from_matrix(b, out);

    matrix_disalloc(b);
    matrix_disalloc(v);
    vector_disalloc(lambda);

    return;
} /* mat_sym2_max_zero_sq */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the absolute value of a symmetric 2*2 matrix. 

    Computes a principal axis transformation and take the absolute value of 
    the eigenvalues. The resulting matrix is calculated from the new 
    eigenvalues and the eigenvectors of the input matrix. 
 */
void mat_sym2_abs(double* in,  /**< input matrix, unchanged */
		  double* out) /**< result, changed         */
{
    double** b;      /* storage for matrix (for PA method) */
    double** v;      /* the eigenvectors */
    double*  lambda; /* the eigenvalues  */

    double delta = 10e-10;
    double eps   = 10e-06;

    matrix_alloc(&b, 4, 4);
    matrix_alloc(&v, 4, 4);
    vector_alloc(&lambda, 4);

    mat_sym2_to_matrix(in, b);

    PA_trans(b, 2, delta, eps, lambda, v);

    lambda[1] = fabs( lambda[1] );
    lambda[2] = fabs( lambda[2] );
 
    PA_backtrans(lambda, v, 2, b);

    mat_sym2_from_matrix(b, out);

    matrix_disalloc(b);
    matrix_disalloc(v);
    vector_disalloc(lambda);

    return;
} /* mat_sym2_abs */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the positive part of a symmetric 2*2 matrix.

    Computes a principal axis transformation and set all negative 
    eigenvalues to zero. The resulting matrix is calculated from the new 
    eigenvalues and the eigenvectors of the input matrix. 
 */
void mat_sym2_pos_part(double* in,  /**< input matrix, unchanged */
		       double* out) /**< result, changed         */
{
    double** b;      /* storage for matrix (for PA method) */
    double** v;      /* the eigenvectors */
    double*  lambda; /* the eigenvalues  */

    double delta = 10e-10;
    double eps   = 10e-06;

    matrix_alloc(&b, 4, 4);
    matrix_alloc(&v, 4, 4);
    vector_alloc(&lambda, 4);

    mat_sym2_to_matrix(in, b);

    PA_trans(b, 2, delta, eps, lambda, v);

    if( lambda[1] < 0.0 )
    {
	lambda[1] = 0.0;
    }
    if( lambda[2] < 0.0 )
    {
	lambda[2] = 0.0;
    }
 
    PA_backtrans(lambda, v, 2, b);

    mat_sym2_from_matrix(b, out);

    matrix_disalloc(b);
    matrix_disalloc(v);
    vector_disalloc(lambda);

    return;
} /* mat_sym2_pos_part */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the maximum of two symmetric 2*2 matrices. 

    The maximum is defined via the formula 
    \f$ \mbox{max}(A, B) := \frac{1}{2} (A + B) + \frac{1}{2} |A - B|\f$ .
 */
void mat_sym2_max(double* in1, /**< first input matrix, unchanged  */
		  double* in2, /**< second input matrix, unchanged */
		  double* out) /**< result, changed                */
{
    double *sum;
    double *diff;
    double *diffabs;

    vector_alloc(&sum, 3);
    vector_alloc(&diff, 3);
    vector_alloc(&diffabs, 3);

    mat_sym2_add_get(in1, in2, sum);
    mat_sym2_sub_get(in1, in2, diff);
    mat_sym2_abs(diff, diffabs);

    mat_sym2_mean_get(sum, diffabs, out);

    vector_disalloc(sum);
    vector_disalloc(diff);
    vector_disalloc(diffabs);

    return;
} /* mat_sym2_max */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the minimum of two symmetric 2*2 matrices. 

    The maximum is defined via the formula 
    \f$ \mbox{min}(A, B) := \frac{1}{2} (A + B) - \frac{1}{2} |A - B|\f$ .
 */
void mat_sym2_min(double* in1, /**< first input matrix, unchanged  */
		  double* in2, /**< second input matrix, unchanged */
		  double* out) /**< result, changed                */
{
    long l;

    double *sum;
    double *diff;
    double *diffabs;

    vector_alloc(&sum, 3);
    vector_alloc(&diff, 3);
    vector_alloc(&diffabs, 3);

    mat_sym2_add_get(in1, in2, sum);
    mat_sym2_sub_get(in1, in2, diff);
    mat_sym2_abs(diff, diffabs);

    for(l=0; l<3; l++)
    {
	out[l] = (sum[l] - diffabs[l]) * 0.5;
    }

    vector_disalloc(sum);
    vector_disalloc(diff);
    vector_disalloc(diffabs);

    return;
} /* mat_sym2_min */

/*--------------------------------------------------------------------------*/

/** \brief Applies a given function f to a symmetric 2*2 matrix.

    Computes a principal axis transformation and applies f to all 
    eigenvalues. Then resulting matrix is calculated from the new 
    eigenvalues and the eigenvectors of the input matrix. 
     
 */
void mat_sym2_function(double* in,          /**< input matrix, unchanged */
		       double (*f)(double,double), 
	                                    /**< the function to apply   */
		       double  param,       /**< parameter for f         */
		       double* out)         /**< result, changed         */
{
    double** b;      /* storage for matrix (for PA method) */
    double** v;      /* the eigenvectors */
    double*  lambda; /* the eigenvalues  */

    double delta = 10e-10;
    double eps   = 10e-06;

    matrix_alloc(&b, 4, 4);
    matrix_alloc(&v, 4, 4);
    vector_alloc(&lambda, 4);

    mat_sym2_to_matrix(in, b);

    PA_trans(b, 2, delta, eps, lambda, v);

    lambda[1] = f( lambda[1], param );
    lambda[2] = f( lambda[2], param );

    PA_backtrans(lambda, v, 2, b);

    mat_sym2_from_matrix(b, out);

    matrix_disalloc(b);
    matrix_disalloc(v);
    vector_disalloc(lambda);

    return;
} /* mat_sym2_function */

/*--------------------------------------------------------------------------*/
/** \brief Determines wheter a symmetric 2*2 matrix is positive 
           semidefinite.

    Returns 1 if the matrix is positive semidefinite, 0 otherwise.
 */
int mat_sym2_positive(double *a)
{
    double **matrix;
    double *lambda;
    double **v;

    double delta = 10e-10;
    double eps   = 10e-06;

    double min_lambda;

    vector_alloc(&lambda, 4);
    matrix_alloc(&matrix, 4, 4);
    matrix_alloc(&v, 4, 4);

    mat_sym2_to_matrix(a, matrix);

    PA_trans(matrix, 2, delta, eps, lambda, v);

    min_lambda = lambda[1];

    if(lambda[2] < min_lambda)
    {
	min_lambda = lambda[2];
    }

    vector_disalloc(lambda);
    matrix_disalloc(matrix);
    matrix_disalloc(v);

    if(min_lambda < 0) 
    {
	return 0;
    }

    return 1;
} /* mat_sym2_positive */

/*--------------------------------------------------------------------------*/
/** \brief Determines wheter a symmetric 2*2 matrix is negative 
           semidefinite.

    Returns 1 if the matrix is negative semidefinite, 0 otherwise.
 */
int mat_sym2_negative(double *a)
{
    double **matrix;
    double *lambda;
    double **v;

    double delta = 10e-10;
    double eps   = 10e-06;

    double max_lambda;

    vector_alloc(&lambda, 4);
    matrix_alloc(&matrix, 4, 4);
    matrix_alloc(&v, 4, 4);

    mat_sym2_to_matrix(a, matrix);

    PA_trans(matrix, 2, delta, eps, lambda, v);

    max_lambda = lambda[1];

    if(lambda[2] > max_lambda)
    {
	max_lambda = lambda[2];
    }

    vector_disalloc(lambda);
    matrix_disalloc(matrix);
    matrix_disalloc(v);

    if(max_lambda > 0) 
    {
	return 0;
    }

    return 1;
} /* mat_sym2_negative */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the maximum of three matrices. 
 */
void mat_sym2_max3(double *a,    /**< input matrix, unchanged */
		   double *b,    /**< input matrix 2, unchanged */
		   double *c,    /**< input matrix 3, unchanged */
		   double *out)  /**< result, changed         */
{
    long l;          /* loop variable */

    double *maxab;
    double *maxca;
    double *maxbc;
    double *S;

    double *maxabc;
    double *maxcab;
    double *maxbca;

    double *diff;

    double *lambda;
    double **matrix;
    double **v;
    double min_lambda;

    double delta = 10e-10;
    double eps   = 10e-06;

    /* allocate storage */
    vector_alloc(&maxab, 3);
    vector_alloc(&maxca, 3);
    vector_alloc(&maxbc, 3);
    vector_alloc(&maxabc, 3);
    vector_alloc(&maxcab, 3);
    vector_alloc(&maxbca, 3);
    vector_alloc(&S, 3);
    vector_alloc(&diff, 3);
    vector_alloc(&lambda, 4);
    matrix_alloc(&matrix, 4, 4);
    matrix_alloc(&v, 4, 4);

    /* take pairwise maxima */
    mat_sym2_max(a, b, maxab);
    mat_sym2_max(maxab, c, maxabc);

    mat_sym2_max(c, a, maxca);
    mat_sym2_max(maxca, b, maxcab);

    mat_sym2_max(b, c, maxbc);
    mat_sym2_max(maxbc, a, maxbca);

    /* and the arithmetic mean */
    for(l=0; l<3; l++)
    {
	S[l] = (maxabc[l] + maxcab[l] + maxbca[l]) / 3.0;
    }

    /* calculate eigenvalues of S - a */
    mat_sym2_sub_get(S, a, diff);
    mat_sym2_to_matrix(diff, matrix);
    PA_trans(matrix, 2, delta, eps, lambda, v);
    /* search for minimum */
    min_lambda = lambda[1];
    if( lambda[2] < min_lambda)
    {
	min_lambda = lambda[2];
    }

    /* calculate eigenvalues of S - b */
    mat_sym2_sub_get(S, b, diff);
    mat_sym2_to_matrix(diff, matrix);
    PA_trans(matrix, 2, delta, eps, lambda, v);
    /* search for minimum */
    for(l=1; l<=2; l++)
    {
	if( lambda[l] < min_lambda)
	{
	    min_lambda = lambda[2];
	}
    }    

    /* calculate eigenvalues of S - c */
    mat_sym2_sub_get(S, c, diff);
    mat_sym2_to_matrix(diff, matrix);
    PA_trans(matrix, 2, delta, eps, lambda, v);
    /* search for minimum */
    for(l=1; l<=2; l++)
    {
	if( lambda[l] < min_lambda)
	{
	    min_lambda = lambda[2];
	}
    }    

/*    printf("minimal lambda is %e\n", min_lambda); */

    /* subtract minimal eigenvalue times the unit matrix from S */
    for(l=0; l<2; l++)
    {
	out[l] = S[l] - min_lambda;
    }
    out[2] = S[2];

    vector_disalloc(maxab);
    vector_disalloc(maxca);
    vector_disalloc(maxbc);
    vector_disalloc(maxabc);
    vector_disalloc(maxcab);
    vector_disalloc(maxbca);
    vector_disalloc(S);
    vector_disalloc(diff);
    vector_disalloc(lambda);
    matrix_disalloc(matrix);
    matrix_disalloc(v);

    return;
} /* mat_sym2_max3 */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the minimum of three matrices. 
 */
void mat_sym2_min3(double *a,    /**< input matrix, unchanged */
		   double *b,    /**< input matrix 2, unchanged */
		   double *c,    /**< input matrix 3, unchanged */
		   double *out)  /**< result, changed         */
{
    long l;          /* loop variable */

    double *minab;
    double *minca;
    double *minbc;
    double *S;

    double *minabc;
    double *mincab;
    double *minbca;

    double *diff;

    double *lambda;
    double **matrix;
    double **v;
    double max_lambda;

    double delta = 10e-10;
    double eps   = 10e-06;

    /* allocate storage */
    vector_alloc(&minab, 3);
    vector_alloc(&minca, 3);
    vector_alloc(&minbc, 3);
    vector_alloc(&minabc, 3);
    vector_alloc(&mincab, 3);
    vector_alloc(&minbca, 3);
    vector_alloc(&S, 3);
    vector_alloc(&diff, 3);
    vector_alloc(&lambda, 4);
    matrix_alloc(&matrix, 4, 4);
    matrix_alloc(&v, 4, 4);

    /* take pairwise minima */
    mat_sym2_min(a, b, minab);
    mat_sym2_min(minab, c, minabc);

    mat_sym2_min(c, a, minca);
    mat_sym2_min(minca, b, mincab);

    mat_sym2_min(b, c, minbc);
    mat_sym2_min(minbc, a, minbca);

    /* and the arithmetic mean */
    for(l=0; l<3; l++)
    {
	S[l] = (minabc[l] + mincab[l] + minbca[l]) / 3.0;
    }

    /* calculate eigenvalues of S - a */
    mat_sym2_sub_get(S, a, diff);
    mat_sym2_to_matrix(diff, matrix);
    PA_trans(matrix, 2, delta, eps, lambda, v);
    /* search for minimum */
    max_lambda = lambda[1];
    if( lambda[2] > max_lambda)
    {
	max_lambda = lambda[2];
    }

    /* calculate eigenvalues of S - b */
    mat_sym2_sub_get(S, b, diff);
    mat_sym2_to_matrix(diff, matrix);
    PA_trans(matrix, 2, delta, eps, lambda, v);
    /* search for minimum */
    for(l=1; l<=2; l++) 
    {
	if( lambda[l] > max_lambda)
	{
	    max_lambda = lambda[l];
	}
    }
    
    /* calculate eigenvalues of S - c */
    mat_sym2_sub_get(S, c, diff);
    mat_sym2_to_matrix(diff, matrix);
    PA_trans(matrix, 2, delta, eps, lambda, v);
    /* search for minimum */
    for(l=1; l<=2; l++) 
    {
	if( lambda[l] > max_lambda)
	{
	    max_lambda = lambda[l];
	}
    }

/*    printf("maximal lambda is %e\n", max_lambda); */

    /* add maximal eigenvalue times the unit matrix to S */
    for(l=0; l<2; l++)
    {
	out[l] = S[l] + max_lambda;
    }
    out[2] = S[2];

    vector_disalloc(minab);
    vector_disalloc(minca);
    vector_disalloc(minbc);
    vector_disalloc(minabc);
    vector_disalloc(mincab);
    vector_disalloc(minbca);
    vector_disalloc(S);
    vector_disalloc(diff);
    vector_disalloc(lambda);
    matrix_disalloc(matrix);
    matrix_disalloc(v);

    return;
} /* mat_sym2_min3 */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the trace of a symmetric 2*2 matrix.
 */
double mat_sym2_trace(double *a)
{
    return a[0] + a[1];
} /* mat_sym2_trace */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the minmod function of three matrices. 
 */
void mat_sym2_minmod_trace(double *a,
			   double *b,
			   double *c,
			   double *out)
{
    long l;

    if( (mat_sym2_trace(a) >= 0.0) && 
	(mat_sym2_trace(b) >= 0.0) &&
	(mat_sym2_trace(c) >= 0.0)    )
	/* all matrices positive semi-definite */
    {
	mat_sym2_min3(a, b, c, out);
    } 
    else if( (mat_sym2_trace(a) <= 0.0) && 
	     (mat_sym2_trace(b) <= 0.0) &&
	     (mat_sym2_trace(c) <= 0.0)    )
        /* all matrices negative semi-definite */
    {
	mat_sym2_max3(a, b, c, out);
    } 
    else               
    {
	for(l=0; l<3; l++)
	{
	    out[l] = 0.0;
	}
    }

    return;
} /* mat_sym2_minmod_trace */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the minmod function of three matrices. 
 */
void mat_sym2_minmod(double *a,
		     double *b,
		     double *c,
		     double *out)
{
    long l;

    if( mat_sym2_positive(a) && mat_sym2_positive(b) 
	&& mat_sym2_positive(c) )
	/* all matrices positive semi-definite */
    {
	mat_sym2_min3(a, b, c, out);
    } 
    else if( mat_sym2_negative(a) && mat_sym2_negative(b) 
	     && mat_sym2_negative(c) ) 
        /* all matrices negative semi-definite */
    {
	mat_sym2_max3(a, b, c, out);
    } 
    else               
    {
	for(l=0; l<3; l++)
	{
	    out[l] = 0.0;
	}
    }

    return;
} /* mat_sym2_minmod */

/*--------------------------------------------------------------------------*/

/** \brief Rotates a symmetric 2x2 tensor by a given angle.

    It is assumed that res = U A U^T, that means the corresponding 
    ellipsoid is rotated in the mathematical positive sense. 
 */
void mat_sym2_rot(double *a,   /**< matrix, input            */
		  double phi,  /**< rotation angle phi       */
		  double *res) /**< resulting matrix, output */
{
    double s = sin(phi);
    double c = cos(phi);

    double c_sqr = c*c;
    double cs    = s*c;
    double s_sqr = s*s;
    
    res[0] = c_sqr*a[0] - 2.0*cs*a[2] + s_sqr*a[1];
    res[1] = s_sqr*a[0] + 2.0*cs*a[2] + c_sqr*a[1];
    res[2] = cs*(a[0]-a[1]) + a[2]*(c_sqr-s_sqr);

    return;
} /* mat_sym2_rot */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the absolute value of the difference of two 
           symmetric 2x2 matrices. 

    The absolute value is to be understood in the matrix-valued sense.
 */
void mat_sym2_abs_diff(double *a,   /**< first matrix, input   */
		       double *b,   /**< second matrix, input  */
		       double *res) /**< result, output        */
{
    mat_sym2_sub_get(a, b, res);
    mat_sym2_abs(res, res);

    return;
} /* mat_sym2_abs_diff */

/*--------------------------------------------------------------------------*/

#endif /* __MAT_SYM2_H__ */
