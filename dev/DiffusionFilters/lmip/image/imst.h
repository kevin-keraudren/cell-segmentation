/*--------------------------------------------------------------------------*/
/** \file imst.h 

   \brief Structure tensor for 2D grey value images.

   Maintainer: Stephan Didas <didas@mia.uni-saarland.de>                    */
/*--------------------------------------------------------------------------*/

#ifndef __IMST_H__
#define __IMST_H__

#include <math.h>

#include "image/immem.h"
#include "matrix/mat_sym2.h"
#include "tensor2d/ten2dconv.h"

/*--------------------------------------------------------------------------*/

/** \brief Finite difference scheme used for the approximation of first 
           derivatives.
 */
enum fd_scheme
{
    CENTRAL,   /**< Central differences */
    SOBEL,     /**< Sobel operators     */
    KUMAR      /**< Kumar operators     */
};


/*--------------------------------------------------------------------------*/

/** \brief Calculates the structure tensor for a 2D grey value image.
 */
void image_structure_tensor(
    double **f,       /**< the image, unchanged                 */
    long   nx,        /**< image size in x direction            */
    long   ny,        /**< image size in y direction            */
    long   bx,        /**< boundary size in x direction         */
    long   by,        /**< boundary size in y direction         */
    double hx,        /**< spatial step size in x direction     */
    double hy,        /**< spatial step size in y direction     */
    double sigma,     /**< pre-smoothing or noise scale         */
    enum fd_scheme scheme, /**< derivative approximation scheme */
    double rho,       /**< post-smoothing or integration scale  */
    double prec,      /**< cutoff precision for convolutions    */
    double ***st)     /**< the structure tensor, output         */
{
    long i, j;        /* loop variables */
    double **v;       /* pre-smoothed version of f */

    double dv_dx;     /* derivative approximation */
    double dv_dy;     /* derivative approximation */

    double w1, w2;    /* weights used for derivative approximation */
    double w3, w4;    /* weights used for derivative approximation */

    /* check boundary size */
    if( bx < 1 || by < 1 )
    {
	fprintf(stderr, "ten2d_structure_tensor: boundary size must be >=1"
		"in all directions.\n");
	fprintf(stderr, "ten2d_structure_tensor: bx: %ld, by: %ld\n",
		bx, by);
	exit(1);
    }

    /* allocate memory for v */
    image_alloc(&v, nx, ny, bx, by);

    /* copy f into v */
    image_copy(f, nx, ny, bx, by, v);

    /* pre-smoothing */
    if(sigma > 0.0)
    {
	image_conv_gauss(v, nx, ny, bx, by, hx, hy, MIRRORING, sigma, prec);
    }

    /* ---- calculate the entries ---- */
    /* mirror boundaries of v */
    image_mirror_bd(v, nx, ny, bx, by);

    switch(scheme)
    {
	case CENTRAL:
	    /* initialisation of weights for central differences */
	    w1 = 1.0 / (2.0*hx);
	    w2 = 1.0 / (2.0*hy);
	    
	    for(i=bx; i<nx+bx; i++)
	    {
		for(j=by; j<ny+by; j++)
		{
		    /* compute the derivatives using central differences */
		    dv_dx = w1 * (  v[i+1][j] - v[i-1][j] );
		    dv_dy = w2 * (  v[i][j+1] - v[i][j-1] );
		    
		    /* calculate entries of the tensor */
		    st[i][j][0] = dv_dx * dv_dx;
		    st[i][j][1] = dv_dy * dv_dy;
		    st[i][j][2] = dv_dx * dv_dy;
		}
	    }
	    break;

	case SOBEL:
	    /* initialisation of weights for Sobel operators */
	    w1 = 1.0 / (8.0*hx);
	    w2 = 1.0 / (4.0*hx);
	    w3 = 1.0 / (8.0*hy);
	    w4 = 1.0 / (4.0*hy);
	    
	    for(i=bx; i<nx+bx; i++)
	    {
		for(j=by; j<ny+by; j++)
		{
		    /* compute the derivatives using Sobel operators */
		    dv_dx = w1 * (  v[i+1][j+1] - v[i-1][j+1]
				    + v[i+1][j-1] - v[i-1][j-1] )
			+   w2 * (  v[i+1][j  ] - v[i-1][j  ] );
		    
		    dv_dy = w3 * (  v[i+1][j+1] - v[i+1][j-1]
				    + v[i-1][j+1] - v[i-1][j-1] )
			+   w4 * (  v[i  ][j+1] - v[i  ][j-1] );
		    
		    /* calculate entries of the tensor */
		    st[i][j][0] = dv_dx * dv_dx;
		    st[i][j][1] = dv_dy * dv_dy;
		    st[i][j][2] = dv_dx * dv_dy;
		}
	    }
	    break;

	case KUMAR:
	    /* initialisation of weights for Kumar operators */
	    w1 = 1.0 / (12.0*hx);
	    w2 = 1.0 / (3.0*hx);
	    w3 = 1.0 / (12.0*hy);
	    w4 = 1.0 / (3.0*hy);
	    
	    for(i=bx; i<nx+bx; i++)
	    {
		for(j=by; j<ny+by; j++)
		{
		    /* compute the derivatives using Kumar operators */
		    dv_dx = w1 * (  v[i+1][j+1] - v[i-1][j+1]
				    + v[i+1][j-1] - v[i-1][j-1] )
			+   w2 * (  v[i+1][j  ] - v[i-1][j  ] );
		    
		    dv_dy = w3 * (  v[i+1][j+1] - v[i+1][j-1]
				    + v[i-1][j+1] - v[i-1][j-1] )
			+   w4 * (  v[i  ][j+1] - v[i  ][j-1] );
		    
		    /* calculate entries of the tensor */
		    st[i][j][0] = dv_dx * dv_dx;
		    st[i][j][1] = dv_dy * dv_dy;
		    st[i][j][2] = dv_dx * dv_dy;
		}
	    }
	    break;

	default:
	    fprintf(stderr, "image_structure_tensor: finite difference "
		    "scheme %d not known.\n", scheme);
	    exit(1);
    }

    /* post-smoothing */
    if(rho > 0.0)
    {
	ten2d_conv_gauss(st, nx, ny, bx, by, 3, hx, hy, DIRICHLET, 
			 rho, prec);
    }

    image_disalloc(v);

    return;
} /* image_structure_tensor */

/*--------------------------------------------------------------------------*/

#endif /* __IMST_H__ */

