/*--------------------------------------------------------------------------*/
/** \file immem.h

   \brief Memory allocation methods for 2D double images.

   Maintainer: Stephan Didas <didas@mia.uni-saarland.de>                    */
/*--------------------------------------------------------------------------*/

#ifndef __IMMEM_H__
#define __IMMEM_H__

#include <stdlib.h>

#include "base/memory.h"

/*--------------------------------------------------------------------------*/

/** \brief Allocates storage for an image of size nx*ny with boundary size 
    bx and by.

    Uses matrix_alloc to allocate storage for a matrix of size
    (nx + 2*bx) * (ny + 2*by) * sizeof(double).

    Exits the program with value 1 if there is not enough storage available.
 */
void image_alloc(double ***image,  /**< the pointer to the image is 
				        stored here                  */
                 long   nx,        /**< image size in x direction    */
                 long   ny,        /**< image size in y direction    */
		 long   bx,        /**< boundary size in x direction */
		 long   by)        /**< boundary size in y direction */
{
    matrix_alloc(image, nx + 2*bx, ny + 2*by);
    
    return;
} /* image_alloc */

/*--------------------------------------------------------------------------*/

/** \brief Disallocates storage of an image. 

    Uses matrix_disalloc to free the storage.
 */
void image_disalloc(double **image)  /**< pointer to the image         */
{
    matrix_disalloc(image);

    return;
} /* image_disalloc */

/*--------------------------------------------------------------------------*/

/** \brief Copies one image into another with same boundary size.
 */
void image_copy(double **source, /**< source image                 */
		long   nx,       /**< image size in x direction    */
		long   ny,       /**< image size in y direction    */
		long   bx,       /**< boundary size in x direction */
		long   by,       /**< boundary size in y direction */
		double **dest)   /**< destination image            */
{
    long i, j;

    for(i=bx; i<nx+bx; i++)
    {
	for(j=by; j<ny+by; j++)
	{
	    dest[i][j] = source[i][j];
	}
    }

    return;
} /* image_copy */

/*--------------------------------------------------------------------------*/

/** \brief Copies one image into another and changes the boundary size.
 */
void image_copy_bd(double **source,  /**< source image                   */
		   long   nx,        /**< image size in x direction      */
		   long   ny,        /**< image size in y direction      */
		   long   bx_source, /**< boundary size of source signal */
		                     /*   in x direction                 */
		   long   by_source, /**< boundary size of source signal */
                                     /*   in y direction                 */
		   long   bx_dest,   /**< boundary size of destination   */
				     /*   image in x direction           */
		   long   by_dest,   /**< boundary size of destination   */
				     /*   image in y direction           */
		   double **dest)    /**< destination image              */
{
    long i, j;

    for(i=0; i<nx; i++)
    {
	for(j=0; j<ny; j++)
	{
	    dest[i+bx_dest][j+by_dest] = source[i+bx_source][j+by_source];
	}
    }

    return;
} /* image_copy_bd */

/*--------------------------------------------------------------------------*/

/** \brief Mirrors the boundary of a 2D image.
 */
void image_mirror_bd(
    double **u,       /**< the image                    */
    long   nx,        /**< image size in x direction    */
    long   ny,        /**< image size in y direction    */
    long   bx,        /**< boundary size in x direction */
    long   by)        /**< boundary size in y direction */
{
    long i, j;

    for(i=bx; i<nx+bx; i++)
    {
	for(j=0; j<by; j++)
	{
	    u[i][by-j-1 ] = u[i][by+j     ];
	    u[i][ny+by+j] = u[i][ny+by-j-1];
	}
    }

    for(i=0; i<bx; i++)
    {
	for(j=0; j<ny+2*by; j++)
	{
	    u[bx-i-1 ][j] = u[bx+i     ][j];
	    u[nx+by+i][j] = u[nx+bx-i-1][j];
	}
    }

    return;
} /* image_mirror_bd */

/*--------------------------------------------------------------------------*/

/** \brief Periodic left (and up) shift of an image, the result is stored 
           into another image.
 */
void image_plshift_get(
    double **in,    /**< input image, unchanged       */
    long   nx,      /**< image size in x direction    */
    long   ny,      /**< image size in y direction    */
    long   bx,      /**< boundary size in x direction */
    long   by,      /**< boundary size in y direction */
    long   shiftx,  /**< displacement in x direction  */
    long   shifty,  /**< displacement in y direction  */
    double **out)   /**< resulting image, output      */
{
    long i, j;
    long distx = shiftx % nx;
    long disty = shifty % ny;

    if( distx < 0 )
    {
	distx += nx;
    }
    
    if( disty < 0 )
    {
	disty += ny;
    }

    for(i=bx+distx; i<nx+bx; i++)
    {
	for(j=by+disty; j<ny+by; j++)
	{
	    out[i-distx][j-disty] = in[i][j];
	}
	
	for(j=by; j<by+disty; j++)
	{
	    out[i-distx][j-disty+ny] = in[i][j];
	}
    }

    for(i=bx; i<bx+distx; i++)
    {
	for(j=by+disty; j<ny+by; j++)
	{
	    out[i-distx+nx][j-disty] = in[i][j];
	}
	
	for(j=by; j<by+disty; j++)
	{
	    out[i-distx+nx][j-disty+ny] = in[i][j];
	}
    }

    return;
} /* image_plshift_get */

/*--------------------------------------------------------------------------*/

/** \brief Periodic left (and up) shift of an image.
 */
void image_plshift(
    double **in,    /**< input image, unchanged       */
    long   nx,      /**< image size in x direction    */
    long   ny,      /**< image size in y direction    */
    long   bx,      /**< boundary size in x direction */
    long   by,      /**< boundary size in y direction */
    long   shiftx,  /**< displacement in x direction  */
    long   shifty)  /**< displacement in y direction  */
{
    double **result;

    image_alloc(&result, nx, ny, bx, by);

    image_plshift_get(in, nx, ny, bx, by, shiftx, shifty, result);

    image_copy(result, nx, ny, bx, by, in);

    image_disalloc(result);

    return;
} /* image_plshift */

/*--------------------------------------------------------------------------*/

#endif /* __IMMEM_H__ */

