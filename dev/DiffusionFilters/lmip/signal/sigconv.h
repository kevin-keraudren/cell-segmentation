/*--------------------------------------------------------------------------*/
/** \file sigconv.h 

   \brief Convolution of 1D double signals. 

   Maintainer: Stephan Didas <didas@mia.uni-saarland.de>                    */
/*--------------------------------------------------------------------------*/

#ifndef __SIGCONV_H__
#define __SIGCONV_H__

#include <math.h>

#include "signal/sigmem.h"

/*--------------------------------------------------------------------------*/

/** \brief The boundary conditions used for the convolution.
 */
enum conv_bc
{ DIRICHLET,  /**< Implicitly extends the signal with zeros in both 
		   directions. */
  MIRRORING,  /**< Mirrors the pixels at the boundaries, e.g. at the left
                   boundary (..., f_2, f_1, f_0, f_1, f_2, ...) is used.   */
  PERIODIC    /**< Extends the signal periodic, e.g.
	           (... f_{n-1}, f_n, f_0, f_1, ...)                       */
};

/*--------------------------------------------------------------------------*/

/** \brief Calculates the convolution of a signal f with a kernel k 
           with periodic boundary conditions.
 */
void signal_conv_per_get_simple(
    double *f,       /**< the signal, input                */
    long   nx,       /**< signal length                    */
    long   bx,       /**< boundary size of the signal      */
    double *kernel,  /**< convolution kernel               */
    long   length,   /**< length of the convolution kernel */
    long   offset,   /**< offset of the kernel index       */
    double *result)  /**< the resulting signal, output     */
{
    long   i, j;     /* loop variables     */
    long   findex;   /* index for signal f */
    double help;     /* variable for sum   */
    
    /* calculate convolution for all pixels */
    for(i=bx; i<nx+bx; i++)
    {
	help   = 0.0;             /* initialise sum */
	findex = i+offset;        /* initialise index for f */
        
        /* periodic boundary conditions */
	while( findex >= nx+bx )  
	{
	    findex = findex - nx;
	}

	/* calculate (f*k)_i */
	for(j=0; j<length; j++)
	{
	    help = help + f[findex]*kernel[j];
	    findex = findex - 1;
	    /* periodic boundary conditions */
	    if(findex < bx)     
	    {
		findex = findex + nx; 
	    } 
	} 
	result[i] = help; 
    }

    return;
} /* signal_conv_per_get_simple */

/*--------------------------------------------------------------------------*/

/** \brief Convolution of an 1d signal with a symmetric kernel.

    \todo Provide an optimised version which avoids allocation of help 
          signal.
 */
void signal_conv_symm(double *f,       /**< initial image, unchanged */
		      long   nx,       /**< image size in x direction */
		      long   bx,       /**< boundary size in x direction */
		      enum conv_bc bc, /**< boundary conditions */
		      double *kernel,  /**< kernel values */
		      long   size,     /**< size of the kernel */
		      double *u)       /**< output image */
{
    long i, k;

    double *help;

    if(size > nx) 
    {
	printf("kernel too large for convolution! \n");
	exit(1);
    }

    signal_alloc(&help, nx, size);

    /* copy f into help */
    signal_copy_bd(f, nx, bx, size, help);

    /* initialise boundary */
    switch(bc)
    {
	case DIRICHLET:
	    for(i=1; i<=size; i++)
	    {
		help[size-i]      = 0.0;
		help[nx+size-1+i] = 0.0;
	    }
	    break;

	case MIRRORING:
	    for(i=1; i<=size; i++)
	    {
		help[size-i]      = help[size+i-1];
		help[nx+size-1+i] = help[nx+size-i];
	    }
	    break;

	case PERIODIC:
	    for(i=1; i<=size; i++)
	    {
		help[size-i]      = help[nx+size-i];
		help[nx+size-1+i] = help[size+i-1];
	    }
	    break;

	default:
	    printf("signal_conv_symm: boundary conditions not supported.\n");
    }


    /* convolution (use symmetry of the kernel) */
    for(i=bx; i<nx+bx; i++)
    {
        u[i] = kernel[0] * help[i-bx+size];
        
        for(k = 1; k<=size; k++)
        {
            u[i] += kernel[k] * (help[i-bx+size-k] + help[i-bx+size+k]);
        }
    }

    signal_disalloc(help);

    return;
} /* signal_conv_symm */

/*--------------------------------------------------------------------------*/

/** \brief Creates a Gaussian signal with given standard deviation.

    Since the Gaussian kernel is symmetric, we only store one half of the 
    values. 
    The kernel is evaluated at equidistant points with spatial step size 
    \f$h_x\f$ and truncated at the point prec * sigma / hx.
    The function allocates a signal \f$ f \f$ with 
    \f$n := \f$(long)(prec * sigma / hx) + 2 
    components and fills in the values of the Gaussian function.
    The values are given as
    \f[
       \tilde{f}(k) := 
       \frac{1}{\sqrt{2 \pi \sigma^2}}
       \exp \left( - \frac{k^2 h^2}{2 \sigma^2} \right)
       \quad \mbox{for all }k \in \{0, \dots, n-1\}
    \f]
    Then
    \f[
       l := \tilde{f}(0) + 2 \sum_{k=1}^{n-1} 2 \tilde{f}(k)
    \f]
    is the 1-norm of the whole truncated Gaussian.
    The output vector is then the normalised version of \f$\tilde{f}\f$,
    i.e. \f$ f := \frac{1}{l} \tilde{f} \f$.
 */
void signal_create_gauss(double **gauss, /**< Gaussian signal    */
			 double hx,      /**< spatial step size  */
			 double sigma,   /**< standard deviation */
			 double prec,    /**< cut-off precision  */
			 long   *size)   /**< size of the kernel */
{
    long   i;

    double argden = hx*hx / (2.0 * sigma * sigma);  
    /* denominator in the argument of the exp function */

    double den    = 1.0 / (sqrt(2 * M_PI) * sigma);
    /* factor in front of the exp function */

    double norm   = 0.0;
    /* to normalise the vector */

    if(sigma <= 0)
    {
	fprintf(stderr, "signal_create_gauss: sigma is %f\n", sigma);
	fprintf(stderr, "signal_create_gauss: it must be larger than "
		"zero.\n");
	exit(1);
    }

    *size = (long)(prec * sigma / hx) + 1;
    /* size of the vector minus 1 */

    /* allocate memory */
    signal_alloc(gauss, *size+1, 0);

    /* initialise Gaussian */
    for(i=0; i <= *size; i++) 
    {
        (*gauss)[i] = den * exp( - i*i * argden );
    }

    /* determine 1-norm */
    norm = (*gauss)[0];
    for(i=1; i <= *size; i++) 
    {
        norm = norm + 2.0 * (*gauss)[i];
    }

    /* and normalise the truncated Gaussian */
    norm = 1.0 / norm;
    for(i=0; i <= *size; i++) 
    {
        (*gauss)[i] = (*gauss)[i] * norm;
    }

    return;
}

/*--------------------------------------------------------------------------*/

/** \brief Convolution of an 1d signal with a Gaussian kernel.

    Uses a Gaussian kernel with given standard deviation. 
 */
void signal_conv_gauss(double *f,       /**< initial image, unchanged */
		       long   nx,       /**< signal length */
		       long   bx,       /**< boundary size */
		       double hx,       /**< spatial step size */
		       double sigma,    /**< standard deviation of Gaussian */
		       double prec,     /**< precision for Gaussian */
		       enum conv_bc bc, /**< boundary conditions */
		       double *u)       /**< resulting image */
{
    double *gauss;  /* vector with Gaussian kernel */
    long size;      /* size of the signal          */

    signal_create_gauss(&gauss, hx, sigma, prec, &size);

    signal_conv_symm(f, nx, bx, bc, gauss, size, u); 

    signal_disalloc(gauss);

    return;
} /* signal_conv_gauss */

/*--------------------------------------------------------------------------*/

#endif /* __SIGCONV_H__ */
