/*--------------------------------------------------------------------------*/
/** \file ten2dmem.h
                                                                              
   \brief Memory allocation methods for 2D tensor fields.   

   Maintainer: Stephan Didas <didas@mia.uni-saarland.de>                    */
/*--------------------------------------------------------------------------*/

#ifndef __TEN2DMEM_H__
#define __TEN2DMEM_H__

#include <stdlib.h>

#include "base/memory.h"

/*--------------------------------------------------------------------------*/

/** \brief Allocates storage for a 2D tensor field of size nx*ny 
    with boundary size bx and by.

    The tensors have dim entries. The allocated memory has size 
    (nx + 2*bx)*(ny + 2*by) * dim * sizeof(double)
    and is allocated as one linear block. 
    
    Exits the program with value 1 if there is not enough storage 
    available.
 */
void ten2d_alloc(double ****tensor, 
                 /**< stores the pointer to the tensor field        */
                 long   nx,   /**< field size in x direction        */
                 long   ny,   /**< field size in y direction        */
		 long   bx,   /**< boundary size in x direction     */
		 long   by,   /**< boundary size in y direction     */
		 long   dim)  /**< number of entries for one tensor */
{
    array3_alloc(tensor, nx + 2*bx, ny + 2*by, dim);
    
    return;
} /* ten2d_alloc */

/*--------------------------------------------------------------------------*/

/** \brief Disallocates storage of a 2D tensor field.
 */
void ten2d_disalloc(double ***tensor) /**< the tensor field         */
{
    array3_disalloc(tensor);

    return;
} /* ten2d_disalloc */

/*--------------------------------------------------------------------------*/

/** \brief Copies one tensor field into another one with same boundary size.
 */
void ten2d_copy(double ***source,  /**< source tensor field             */
		long   nx,         /**< field size in x direction       */
		long   ny,         /**< field size in y direction       */
		long   bx,         /**< boundary size in x direction    */
		long   by,         /**< boundary size in y direction    */
                long   dim,        /**< number of entries in one tensor */
		double ***dest)    /**< destination tensor field        */
{
    long i, j, l; /* loop variables */

    for(i=bx; i<nx+bx; i++)
    {
	for(j=by; j<ny+by; j++)
	{
	    for(l=0; l<dim; l++)
	    {
		dest[i][j][l] = source[i][j][l];
	    }
	}
    }

    return;
} /* ten2d_copy */

/*--------------------------------------------------------------------------*/

/** \brief Copies one 2D tensor field into another one and changes the 
           boundary size. 
 */
void ten2d_copy_bd(double ****source, /**< source tensor field              */
		   long   nx,         /**< field size in x direction        */
		   long   ny,         /**< field size in y direction        */
		   long   bx_source,  /**< boundary size of source field    */
		                      /*   in x direction                   */
		   long   by_source,  /**< boundary size of source field    */
                                      /*   in y direction                   */
		   long   dim,        /**< number of entries for one tensor */
		   long   bx_dest,    /**< boundary size of destination     */
				      /*   field in x direction             */
		   long   by_dest,    /**< boundary size of destination     */
				      /*   field in y direction             */
		   double ****dest)   /**< destination tensor field         */
{
    long i, j, l;

    for(i=0; i<nx; i++)
    {
	for(j=0; j<ny; j++)
	{
	    for(l=0; l<dim; l++)
	    {
		dest[i+bx_dest][j+by_dest][l] 
		    = source[i+bx_source][j+by_source][l];
	    }
	}
    }

    return;
} /* ten2d_copy_bd */

/*--------------------------------------------------------------------------*/

/** \brief Fills the boundary of a tensor field by mirroring the entries.
 */
void ten2d_mirror_bd(
    double ***u,   /**< the tensor field                 */
    long   nx,     /**< field size in x direction        */
    long   ny,     /**< field size in y direction        */
    long   bx,     /**< boundary size in x direction     */
    long   by,     /**< boundary size in y direction     */
    long   dim)    /**< number of entries for one tensor */
{
    long i, j, l;

    /* mirror in y direction */
    for(i=bx; i<nx+bx; i++)
    {
	for(j=0; j<by; j++)
	{
	    for(l=0; l<dim; l++)
	    {
		u[i][by-j-1 ][l] = u[i][by+j     ][l];
		u[i][ny+by+j][l] = u[i][ny+by-j-1][l];
	    }
	}
    }

    /* mirror in x direction */
    for(i=0; i<bx; i++)
    {
	for(j=0; j<ny+2*by; j++)
	{
	    for(l=0; l<dim; l++)
	    {
		u[bx-i-1 ][j][l] = u[bx+i     ][j][l];
		u[nx+bx+i][j][l] = u[nx+bx-i-1][j][l];
	    }
	}
    }

    return;
} /* ten2d_mirror_bd */

/*--------------------------------------------------------------------------*/

/** \brief Periodic left (and up) shift of a 2D tensor field, the result 
           is stored into another field.

    ATTENTION: untested
 */
void ten2d_plshift_get(
    double ***in,   /**< input field, unchanged       */
    long   nx,      /**< image size in x direction    */
    long   ny,      /**< image size in y direction    */
    long   bx,      /**< boundary size in x direction */
    long   by,      /**< boundary size in y direction */
    long   dim,     /**< number of entries per tensor */
    long   shiftx,  /**< displacement in x direction  */
    long   shifty,  /**< displacement in y direction  */
    double ***out)  /**< resulting field, output      */
{
    long i, j, l;   /* loop variables */
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
            for(l=0; l<dim; l++)
            {
                out[i-distx][j-disty][l] = in[i][j][l];
            }
        }

        for(j=by; j<by+disty; j++)
        {
            for(l=0; l<dim; l++)
            {
                out[i-distx][j-disty+ny][l] = in[i][j][l];
            }
        }
    }

    for(i=bx; i<bx+distx; i++)
    {
        for(j=by+disty; j<ny+by; j++)
        {
            for(l=0; l<dim; l++)
            {
                out[i-distx+nx][j-disty][l] = in[i][j][l];
            }
        }

        for(j=by; j<by+disty; j++)
        {
            for(l=0; l<dim; l++)
            {
                out[i-distx+nx][j-disty+ny][l] = in[i][j][l];
            }
        }
    }

    return;
} /* ten2d_plshift_get */

/*--------------------------------------------------------------------------*/

/** \brief Periodic left (and up) shift of a 2D tensor field.

    ATTENTION: untested
 */
void ten2d_plshift(
    double ***in,   /**< input image, unchanged       */
    long   nx,      /**< image size in x direction    */
    long   ny,      /**< image size in y direction    */
    long   bx,      /**< boundary size in x direction */
    long   by,      /**< boundary size in y direction */
    long   dim,     /**< number of entries per tensor */
    long   shiftx,  /**< displacement in x direction  */
    long   shifty)  /**< displacement in y direction  */
{
    double ***result;

    ten2d_alloc(&result, nx, ny, bx, by, dim);

    ten2d_plshift_get(in, nx, ny, bx, by, dim, shiftx, shifty, result);

    ten2d_copy(result, nx, ny, bx, by, dim, in);

    ten2d_disalloc(result);

    return;
} /* ten2d_plshift */

/*--------------------------------------------------------------------------*/

#endif /* __TEN2DMEM_H__ */

