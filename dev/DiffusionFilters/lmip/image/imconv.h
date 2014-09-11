/*--------------------------------------------------------------------------*/
/** \file imconv.h 

   \brief Convolution of 2D double images. 

   Maintainer: Stephan Didas <didas@mia.uni-saarland.de>                    */
/*--------------------------------------------------------------------------*/

#ifndef __IMCONV_H__
#define __IMCONV_H__

#include <math.h>

#include "signal/sigmem.h"
#include "signal/sigconv.h"

/*--------------------------------------------------------------------------*/

/** \brief Computes the convolution of all image columns with a given 
           symmetric 1D kernel.

    Uses the function signal_conv_symm for convolution.
 */
void image_conv_symm_cols(
    double  **f,      /**< input image, unchanged       */
    long    nx,       /**< image size in x direction    */
    long    ny,       /**< image size in x direction    */
    long    bx,       /**< boundary size in x direction */
    long    by,       /**< boundary size in y direction */
    enum conv_bc bc,  /**< boundary conditions          */
    double  *kernel,  /**< convolution kernel           */
    long    size,     /**< size of the kernel           */
    double  **u)      /**< output image                 */
{
    long i;  /* loop variable */

    for(i=bx; i<nx+bx; i++)
    {
	signal_conv_symm(f[i], ny, by, bc, kernel, size, u[i]);
    }

    return;
} /* image_conv_symm_cols */

/*--------------------------------------------------------------------------*/

/** \brief Computes the convolution of all image rows with a given symmetric
           1D kernel.

    Uses the function signal_conv_symm for convolution.
 */
void image_conv_symm_rows(
    double  **f,      /**< input image, unchanged       */
    long    nx,       /**< image size in x direction    */
    long    ny,       /**< image size in x direction    */
    long    bx,       /**< boundary size in x direction */
    long    by,       /**< boundary size in y direction */
    enum conv_bc bc,  /**< boundary conditions          */
    double  *kernel,  /**< convolution kernel           */
    long    size,     /**< size of the kernel           */
    double  **u)      /**< output image                 */
{
    long i, j; /* loop variables */
    double *row_in; /* time saver */
    double *row_out;

    signal_alloc(&row_in, nx, bx);
    signal_alloc(&row_out, nx, bx);

    for(j=by; j<ny+by; j++)
    {
	for(i=bx; i<nx+bx; i++)
	{
	    row_in[i] = f[i][j];
	}

	signal_conv_symm(row_in, nx, bx, bc, kernel, size, row_out);

	for(i=bx; i<nx+bx; i++)
	{
	    u[i][j] = row_out[i];
	}
    }

    signal_disalloc(row_in);
    signal_disalloc(row_out);

    return;
} /* image_conv_symm_rows */

/*--------------------------------------------------------------------------*/

/** \brief Computes the convolution of an image with a Gaussian kernel.
 */
void image_conv_gauss(
    double **f,   /**< in- and output image             */
    long   nx,    /**< image size in x direction        */
    long   ny,    /**< image size in y direction        */
    long   bx,    /**< boundary size in x direction     */
    long   by,    /**< boundary size in y direction     */
    double hx,    /**< spatial step size in x direction */
    double hy,    /**< spatial step size in y direction */
    enum conv_bc bc,   /**< boundary conditions              */
    double sigma, /**< standard deviation of Gaussian   */
    double prec)  /**< cutoff precision                 */
{
    double *gauss; /* one half of the Gaussian signal     */
    long   size;   /* size of the half of Gaussian signal */
    double **help;

    image_alloc(&help, nx, ny, bx, by);

    /* convolution in x direction */
    signal_create_gauss(&gauss, hx, sigma, prec, &size);
    image_conv_symm_rows(f, nx, ny, bx, by, bc, gauss, size, help);
    signal_disalloc(gauss);

    /* convolution in y direction */
    signal_create_gauss(&gauss, hy, sigma, prec, &size);
    image_conv_symm_cols(help, nx, ny, bx, by, bc, gauss, size, f);
    signal_disalloc(gauss);

    image_disalloc(help);

    return;
} /* image_conv_gauss */

/*--------------------------------------------------------------------------*/

#endif /* __IMCONV_H__ */
