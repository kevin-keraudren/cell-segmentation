/*--------------------------------------------------------------------------*/
/** \file ten2dconv.h 

   \brief Convolution of 2D tensor fields.

   Maintainer: Stephan Didas <didas@mia.uni-saarland.de>                    */
/*--------------------------------------------------------------------------*/

#ifndef __TEN2DCONV_H__
#define __TEN2DCONV_H__

#include <math.h>

#include "signal/sigmem.h"
#include "signal/sigconv.h"
#include "tensor2d/ten2dmem.h"

/*--------------------------------------------------------------------------*/

/** \brief Computes the convolution of all tensor field columns with a 
           given symmetric 1D kernel.

    Uses the function signal_conv_symm for convolution.
 */
void ten2d_conv_symm_cols(
    double  ***f,    /**< input tensor field, unchanged   */
    long    nx,      /**< image size in x direction       */
    long    ny,      /**< image size in x direction       */
    long    bx,      /**< boundary size in x direction    */
    long    by,      /**< boundary size in y direction    */
    long    dim,     /**< number of entries in one tensor */
    enum conv_bc bc, /**< boundary conditions             */
    double  *kernel, /**< convolution kernel              */
    long    size,    /**< size of the kernel              */
    double  ***u)    /**< output tensor field             */
{
    long i, j, k; /* loop variables */
    double *row_in;   /* time saver */
    double *row_out;  /* time saver */

    signal_alloc(&row_in, ny, by);
    signal_alloc(&row_out, ny, by);

    for(k=0; k<dim; k++)
    {
	for(i=bx; i<nx+bx; i++)
	{
	    for(j=by; j<ny+by; j++)
	    {
		row_in[j] = f[i][j][k];
	    }

	    signal_conv_symm(row_in, ny, by, bc, kernel, size, row_out);
	    
	    for(j=by; j<ny+by; j++)
	    {
		u[i][j][k] = row_out[j];
	    }
	}
    }

    signal_disalloc(row_in);
    signal_disalloc(row_out);

    return;
} /* ten2d_conv_symm_cols */

/*--------------------------------------------------------------------------*/

/** \brief Computes the convolution of all tensor field rows with a 
           given symmetric 1D kernel.

    Uses the function signal_conv_symm for convolution.
 */
void ten2d_conv_symm_rows(
    double  ***f,    /**< input image, unchanged          */
    long    nx,      /**< image size in x direction       */
    long    ny,      /**< image size in x direction       */
    long    bx,      /**< boundary size in x direction    */
    long    by,      /**< boundary size in y direction    */
    long    dim,     /**< number of entries in one tensor */
    enum conv_bc bc, /**< boundary conditions             */
    double  *kernel, /**< convolution kernel              */
    long    size,    /**< size of the kernel              */
    double  ***u)    /**< output image                    */
{
    long i, j, k;    /* loop variables */
    double *row_in;  /* time saver */
    double *row_out; /* time saver */

    signal_alloc(&row_in, nx, bx);
    signal_alloc(&row_out, nx, bx);

    for(k=0; k<dim; k++)
    {
	for(j=by; j<ny+by; j++)
	{
	    for(i=bx; i<nx+bx; i++)
	    {
		row_in[i] = f[i][j][k];
	    }

	    signal_conv_symm(row_in, nx, bx, bc, kernel, size, row_out);
	    
	    for(i=bx; i<nx+bx; i++)
	    {
		u[i][j][k] = row_out[i];
	    }
	}
    }

    signal_disalloc(row_in);
    signal_disalloc(row_out);

    return;
} /* ten2d_conv_symm_rows */

/*--------------------------------------------------------------------------*/

/** \brief Computes the convolution of an 2D tensor field with a Gaussian 
           kernel.
 */
void ten2d_conv_gauss(
    double ***f,     /**< in- and output image             */
    long   nx,       /**< image size in x direction        */
    long   ny,       /**< image size in y direction        */
    long   bx,       /**< boundary size in x direction     */
    long   by,       /**< boundary size in y direction     */
    long   dim,      /**< number of entries in one tensor  */
    double hx,       /**< spatial step size in x direction */
    double hy,       /**< spatial step size in y direction */
    enum conv_bc bc, /**< boundary conditions              */
    double sigma,    /**< standard deviation of Gaussian   */
    double prec)     /**< cutoff precision                 */
{
    double *gauss;   /* one half of the Gaussian signal     */
    long   size;     /* size of the half of Gaussian signal */
    double ***help;  /* time saver */

    ten2d_alloc(&help, nx, ny, bx, by, dim);

    /* convolution in x direction */
    signal_create_gauss(&gauss, hx, sigma, prec, &size);
    ten2d_conv_symm_rows(f, nx, ny, bx, by, dim, bc, gauss, size, help);
    signal_disalloc(gauss);

    /* convolution in y direction */
    signal_create_gauss(&gauss, hy, sigma, prec, &size);
    ten2d_conv_symm_cols(help, nx, ny, bx, by, dim, bc, gauss, size, f);
    signal_disalloc(gauss);

    ten2d_disalloc(help);

    return;
} /* image_conv_gauss */

/*--------------------------------------------------------------------------*/

#endif /* __IMCONV_H__ */
