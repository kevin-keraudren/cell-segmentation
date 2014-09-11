/*--------------------------------------------------------------------------*/
/** \file file.h

   \brief Basic methods and definitions for file handling.

   Maintainer: Stephan Didas <didas@mia.uni-saarland.de>                    */
/*--------------------------------------------------------------------------*/

#ifndef __FILE_H__
#define __FILE_H__

/*--------------------------------------------------------------------------*/

/** \brief Types of image files.
 */
enum file_type
{
    PGM,
    PPM,
    PDM
};

/*--------------------------------------------------------------------------*/

/** \brief Representation of values in PGM and PPM files.
 */
enum pbm_storage_type
{
    RAW,   /**< Pure bindary representation of the grey values.          */
    ASCII  /**< Human readable ASCII representation of quantised integer
	      values.                                                    */
};

/*--------------------------------------------------------------------------*/

/** \brief For PDM files: Representation of the data.
 */
enum pdm_storage_type
{
    F,  /**< Human readable ASCII representation of the floating point 
	   numbers.                                                    */
    I,  /**< Human readable ASCII representation of quantised integer
	   values.                                                     */
    B   /**< Quantisation to integer values between 0 and 255 and 
	   storing as unsigned char per pixel/voxel.                   */
};

/*--------------------------------------------------------------------------*/

#endif /* __FILE_H__ */
