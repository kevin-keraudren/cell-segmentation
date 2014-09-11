/*--------------------------------------------------------------------------*/
/** \file matrix.h   

   \brief Basic operations with n*n double matrices. 

   Maintainer: Stephan Didas <didas@mia.uni-saarland.de>                    */
/*--------------------------------------------------------------------------*/

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <math.h>
#include <stdlib.h>

#include "base/memory.h"

/*--------------------------------------------------------------------------*/

/** \brief Calculates the principal axis transformation of a symmetric 
           n*n matrix. 
	   
   Cyclic Jacobi method for determining the eigenvalues and eigenvectors
   of a symmetric matrix.
   Ref.:  H.R. Schwarz: Numerische Mathematik. Teubner, Stuttgart, 1988.
   pp. 243-246.

   The eigenvectors and eigenvalues are ordered such that eigenvalues 
   are nonincreasing. 

   \note The matrix B is stored with indices 1 ... nx and 1 ... ny.
*/
void PA_trans
(double **B,     /**< in: matrix; only lower triangle needed; unchanged */
 long   n,       /**< B is n * n matrix */
 double delta,   /**< machine precision: lowest delta with 1 + delta > 1 */
 double eps,     /**< absolute precision of the calculated eigenvalues */
 double *lambda, /**< out: eigenvalues, ordered with decreasing modulus */
 double **v)     /**< out: orthonormal eigenvectors as columns */
{
    long   i, j, p, q;    /* loop variables */
    double sum;           /* for summing up */
    double theta, t, c;   /* intermediate results */
    double r, s, g, h;    /* intermediate results */
    long   k;             /* for switching */
    double help;          /* for switching */
    double **a;           /* working copy of B */

    /* ---- allocate storage ---- */
    matrix_alloc (&a, n+1, n+1);

    /* ---- initializations ---- */

    /* a := B */
    for (i=1; i<=n; i++)
    {
	for (j=1; j<=n; j++)
	{
	    a[i][j] = B[i][j];
	}
    }

    /* v := unit matrix */
    for (i=1; i<=n; i++)
    {
	for (j=1; j<=n; j++)
	{
	    v[i][j] = 0.0;
	}
	v[i][i] = 1;
    }


    /* ---- loop ---- */

    do {
	/* check whether accuracy is reached */
	sum = 0.0;
	for (i=2; i<=n; i++)
	{
	    for (j=1; j<=i-1; j++)
	    {
		sum = sum + a[i][j] * a[i][j];
	    }
	}
	
	if (sum + sum > eps * eps)
	    /* accuracy not reached yet, new cycle */
	{
	    for (p=1; p<=n-1; p++)
	    {
		for (q=p+1; q<=n; q++)
		{
		    if (fabs(a[q][p]) >= eps * eps)
		    {
			theta = (a[q][q] - a[p][p]) / (2.0 * a[q][p]);
			t = 1.0;
			if (fabs(theta) > delta)
			    t = 1.0 / (theta + theta/fabs(theta) *
				       sqrt (theta * theta + 1.0));
			c = 1.0 / sqrt (1.0 + t * t);
			s = c * t;
			r = s / (1.0 + c);
			a[p][p] = a[p][p] - t * a[q][p];
			a[q][q] = a[q][q] + t * a[q][p];
			a[q][p] = 0.0;
			for (j=1; j<=p-1; j++)
			{
			    g = a[q][j] + r * a[p][j];
			    h = a[p][j] - r * a[q][j];
			    a[p][j] = a[p][j] - s * g;
			    a[q][j] = a[q][j] + s * h;
			}
			for (i=p+1; i<=q-1; i++)
			{
			    g = a[q][i] + r * a[i][p];
			    h = a[i][p] - r * a[q][i];
			    a[i][p] = a[i][p] - s * g;
			    a[q][i] = a[q][i] + s * h;
			}
			for (i=q+1; i<=n; i++)
			{
			    g = a[i][q] + r * a[i][p];
			    h = a[i][p] - r * a[i][q];
			    a[i][p] = a[i][p] - s * g;
			    a[i][q] = a[i][q] + s * h;
			}
			for (i=1; i<=n; i++)
			{
			    g = v[i][q] + r * v[i][p];
			    h = v[i][p] - r * v[i][q];
			    v[i][p] = v[i][p] - s * g;
			    v[i][q] = v[i][q] + s * h;
			}
		    } /* if */
		} /* for */
	    } /* for */
	} /* if */
    } /* do */
    while (sum + sum > eps * eps);
	    
    for (i=1; i<=n; i++)
    {
	lambda[i] = a[i][i];
    }
    
	    
    /* ---- order eigenvalues and eigenvectors ---- */
    
    for (i=1; i<=n-1; i++)
    {
	k = i;
	for (j=i+1; j<=n; j++)
	{
	    /*  if (fabs(lambda[j]) > fabs(lambda[k])) */
	    if ((lambda[j]) > (lambda[k]))
	    {
		k = j;
	    }
	}

	if (k != i)
	{
	    /* switch eigenvalue i and k */
	    help      = lambda[k];
	    lambda[k] = lambda[i];
	    lambda[i] = help;
	    /* switch eigenvector i and k */
	    for (j=1; j<=n; j++)
	    {
		help    = v[j][k];
		v[j][k] = v[j][i];
		v[j][i] = help;
	    }
	} /* if */
    } /* for */
    
    /* ---- disallocate storage ---- */
    matrix_disalloc(a);

    return;
} /* PA_trans */

/*--------------------------------------------------------------------------*/

/** \brief Calculates a symmetric n*n matrix from its eigenvectors and 
           eigenvalues. 

    \note The output matrix b is stored with indices 1 ... nx and 1 ... ny.
 */
void PA_backtrans(double *lambda, /**< input: eigenvalues              */
		  double **v,     /**< input: orthonormal eigenvectors */
                                  /*   as columns                      */
		  long   n,       /**< matrix has size n * n           */
		  double **b)     /**< out: matrix                     */
{
    long   i, j, k;    /* loop variables */
    double sum;        /* for summing up */
    
    /* calculate lower triangle */
    for (i=1; i<=n; i++)
    {
	for (j=1; j<=i; j++)
	{
	    sum = lambda[1] * v[i][1] * v[j][1];
	    for (k=2; k<=n; k++)
	    {
		sum = sum + lambda[k] * v[i][k] * v[j][k];
	    }
	    b[i][j] = sum;
	}
    }
    
    /* upper triangle by mirroring */
    for (i=1; i<=n; i++)
    {
	for (j=i+1; j<=n; j++)
	{
	    b[i][j] = b[j][i];
	}
    }

    return;
} /* PA_backtrans */

/*--------------------------------------------------------------------------*/

#endif /* __MATRIX_H__ */

