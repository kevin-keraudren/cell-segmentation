/*--------------------------------------------------------------------------*/
/** \file sigmem.h

   \brief Memory allocation methods for 1D double signals.

   Maintainer: Stephan Didas <didas@mia.uni-saarland.de>                    */
/*--------------------------------------------------------------------------*/

#ifndef __SIGMEM_H__
#define __SIGMEM_H__

#include <stdlib.h>

#include "base/memory.h"

/*--------------------------------------------------------------------------*/

/** \brief Allocates storage for a signal of size nx with boundary size bx.

    The allocated memory has size nx + 2*bx. 
    Exits the program with value 1 if there is not enough storage available.
 */
void signal_alloc(double **signal,      /**< the pointer to the signal is 
                                             stored here                   */
                  long   nx,            /**< size of the signal            */
		  long   bx)            /**< size of the boundary          */
{
    vector_alloc(signal, nx + 2*bx);
    
    return;
} /* signal_alloc */

/*--------------------------------------------------------------------------*/

/** \brief Disallocates storage of a signal. 
 */
void signal_disalloc(double *signal)    /**< pointer to the signal  */
{
    vector_disalloc(signal);

    return;
} /* signal_disalloc */

/*--------------------------------------------------------------------------*/

/** \brief Copies one signal into another with same boundary size.
 */
void signal_copy(double *source,        /**< source signal       */
		 long   nx,             /**< size of the signals */
		 long   bx,             /**< boundary size       */
		 double *dest)          /**< destination signal  */
{
    long i;

    for(i=bx; i<nx+bx; i++)
    {
	dest[i] = source[i];
    }

    return;
} /* signal_copy */

/*--------------------------------------------------------------------------*/

/** \brief Copies one signal into another and changes the boundary size.
 */
void signal_copy_bd(double *source,     /**< source signal                  */
		    long   nx,          /**< size of the signals            */
		    long   bx_source,   /**< boundary size of source signal */
		    long   bx_dest,     /**< boundary size of destination 
                                             signal                         */
		    double *dest)       /**< destination signal             */
{
    long i;

    for(i=0; i<nx; i++)
    {
	dest[i+bx_dest] = source[i+bx_source];
    }

    return;
}

/*--------------------------------------------------------------------------*/

/** \brief Copies a small signal at a certain position into a larger one.
 */
void signal_copy_into(
    double *source, /**< source signal                      */
    long   nx,      /**< size of the area to copy           */
    long   bx,      /**< boundary size of the source signal */
    double *dest,   /**< destination signal                 */
    long   start,   /**< position in the destination signal
		         to copy source                     */
    long   bx_dest) /**< boundary size (destination signal) */
{
    long i;

    for(i=bx; i<nx+bx; i++)
    {
	dest[i-bx+bx_dest+start] = source[i];
    }

    return;
} /* signal_copy_into */

/*--------------------------------------------------------------------------*/

/** \brief Periodic left shift of a signal, the result is stored in to 
           another signal.

    shift can be an arbitrary integer number, we only use shift modulo nx
    as displacement.
 */
void signal_plshift_get(
    double *in,     /**< signal, input                         */
    long   nx,      /**< size of the signal                    */
    long   bx,      /**< boundary size of the signal           */
    long   shift,   /**< displacement of the shift to the left */
    double *out)    /**< resulting signal, output              */
{
    long i;
    long dist = shift % nx;

    if( dist < 0 ) {
	dist = dist + nx;
    }

    for(i=bx+dist; i<nx+bx; i++)
    {
	out[i-dist] = in[i];
    }

    for(i=bx; i<bx+dist; i++)
    {
	out[i-dist+nx] = in[i];
    }

    return;
} /* signal_plshift_get */

/*--------------------------------------------------------------------------*/

/** \brief Periodic left shift of a signal.
 */
void signal_plshift(
    double *signal,   /**< signal, in- and output                */
    long   nx,        /**< size of the signal                    */
    long   bx,        /**< boundary size                         */
    long   shift)     /**< displacement of the shift to the left */
{
    double *result;

    signal_alloc(&result, nx, bx);

    signal_plshift_get(signal, nx, bx, shift, result);

    signal_copy(result, nx, bx, signal);

    signal_disalloc(result);

    return;
} /* signal_plshift_get */

/*--------------------------------------------------------------------------*/

#endif /* __SIGMEM_H__ */
