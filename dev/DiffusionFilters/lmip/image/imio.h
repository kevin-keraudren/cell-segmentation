/*--------------------------------------------------------------------------*/
/** \file imio.h

   \brief Input and output methods for 2D double images.

   Maintainer: Stephan Didas <didas@mia.uni-saarland.de>                    */
/*--------------------------------------------------------------------------*/

#ifndef __IMIO_H_
#define __IMIO_H_

#include <stdio.h>
#include <string.h>

#include "base/comment.h"
#include "base/file.h"
#include "image/immem.h"

/*--------------------------------------------------------------------------*/

/** \brief Reads a 2D grey value image from a file in raw PGM format.
    
    Basically follows the specification in the Netpbm package
    http://netpbm.sourceforge.net, but does only implement a limited range 
    of this specification.
 */
void image_read_pgm(
    const char *file_name,   /**< name of file to read                    */
    double ***image, /**< the pointer to the image is stored here         */
    long   *nx,      /**< image size in x direction as given in the input 
                          file                                            */
    long   *ny,      /**< image size in y direction as given in the input 
                          file                                            */
    long   bx,       /**< boundary size in x direction                    */
    long   by,       /**< boundary size in y direction                    */
    struct comment **comments) /**< the comment field of the input file   */
{
    FILE* in_file;
    char  row[256];
    long  i, j;
    int   maxval;

    /* open input file */
    in_file = fopen(file_name, "r");
    if(in_file == 0)
    {
	printf("image_read_pgm: Cannot open input file %s.\n", 
	       file_name);
	exit(1);
    }

    /* check file type with magic number */
    fgets(row, 255, in_file);
    if(strncmp(row, "P5", 2))
    {
	printf("image_read_pgm: File format is not raw PGM.\n");
	exit(1);
    }
    
    /* read comments */
    comment_alloc(comments);
    fgets(row, 255, in_file);
    while(row[0] == '#')
    {
	comment_add(*comments, row);
	fgets(row, 255, in_file);
    }

    /* read width (nx) and height (ny) */
    sscanf(row, "%ld %ld", nx, ny);

    /* read maximal grey value */
    fgets(row, 255, in_file);
    sscanf(row, "%d", &maxval);

    /* allocate memory for the image */
    image_alloc(image, *nx, *ny, bx, by);

    /* read image data */
    if(maxval < 256) {
	/* one byte per pixel */
	for(j=by; j<*ny+by; j++)
	{
	    for(i=bx; i<*nx+bx; i++)
	    {
		(*image)[i][j] = (double) getc(in_file);
	    }
	}
    } else {
	/* two bytes per pixel, most significant byte first */
	for(j=by; j<*ny+by; j++)
	{
	    for(i=bx; i<*nx+bx; i++)
	    {
		(*image)[i][j] = (double) getc(in_file) + 
		    (double) getc(in_file)*256;
	    }
	}
    }

    fclose(in_file);
    printf("image_read_pgm: input file %s read.\n", file_name);

    return;
} /* image_read_pgm */

/*--------------------------------------------------------------------------*/

/** \brief Writes a 2D grey value image to a file in raw PGM format.
    
    Basically follows the specification in the Netpbm package
    http://netpbm.sourceforge.net, but does only implement a limited range 
    of this specification.

    This method only supports a maximal grey value of 255 so far.
    The grey values are truncated to 0 and 255.
 */
void image_write_pgm(
    const char *file_name,     /**< name of file to write        */
    double **image,            /**< image data                   */
    long   nx,                 /**< image size in x direction    */
    long   ny,                 /**< image size in y direction    */
    long   bx,                 /**< boundary size in x direction */
    long   by,                 /**< boundary size in y direction */
    struct comment *comments)  /**< comment field                */
{
    FILE* out_file;
    long  i, j;
    unsigned char byte;

    /* open output file */
    out_file = fopen(file_name, "w");
    if(out_file == 0)
    {
	printf("image_write_pgm: Could not open output file %s.\n", 
	       file_name);
	exit(1);
    }

    /* ---- PGM header ---- */
    fprintf(out_file, "P5\n");

    /* print out comments */
    while(comments != NULL)
    {
	if(comments->line != NULL) 
	{
	    fprintf(out_file, "%s", comments->line);
	}
	comments = comments->next;
    }

    /* image size in x and y direction */
    fprintf(out_file, "%ld %ld\n", nx, ny);

    /* maximal grey value */
    fprintf(out_file, "%d\n", 255);

    /* write image data */
    for(j=by; j<ny+by; j++)
    {
	for(i=bx; i<nx+bx; i++)
	{
	    if(image[i][j] < 0) 
	    {
		byte = (unsigned char) 0.0;
	    } 
	    else if(image[i][j] > 255.0)
	    {
		byte = (unsigned char) 255.0;
	    } 
	    else
	    {
		byte = (unsigned char) image[i][j];
	    }

	    fwrite(&byte, sizeof(unsigned char), 1, out_file);
	}
    }

    fclose(out_file);

    return;
} /* image_write_pgm */

/*--------------------------------------------------------------------------*/
/** \brief Writes a 2D grey value image to a file in raw PGM format with
    16 bits per pixel.

    Basically follows the specification in the Netpbm package
    http://netpbm.sourceforge.net, but does only implement a limited range 
    of this specification.

    This method only supports a maximal grey value of 255 so far.
    The grey values are truncated to 0 and 255.
 */
void image_write_pgm_16bpp(
    const char *file_name,     /**< name of file to write        */
    double **image,            /**< image data                   */
    long   nx,                 /**< image size in x direction    */
    long   ny,                 /**< image size in y direction    */
    long   bx,                 /**< boundary size in x direction */
    long   by,                 /**< boundary size in y direction */
    struct comment *comments)  /**< comment field                */
{
    FILE* out_file;
    long  i, j;
    long  value;
    unsigned char byte;

    /* open output file */
    out_file = fopen(file_name, "w");
    if(out_file == 0)
    {
	printf("image_write_pgm: Could not open output file %s.\n", 
	       file_name);
	exit(1);
    }

    /* ---- PGM header ---- */
    fprintf(out_file, "P5\n");

    /* print out comments */
    while(comments != NULL)
    {
	if(comments->line != NULL) 
	{
	    fprintf(out_file, "%s", comments->line);
	}
	comments = comments->next;
    }

    /* image size in x and y direction */
    fprintf(out_file, "%ld %ld\n", nx, ny);

    /* maximal grey value */
    fprintf(out_file, "%d\n", 65535);

    /* write image data */
    for(j=by; j<ny+by; j++)
    {
	for(i=bx; i<nx+bx; i++)
	{
	    if(image[i][j] < 0) 
	    {
		value = 0;
	    } 
	    else if(image[i][j] > 65535.0)
	    {
		value = 65535;
	    } 
	    else
	    {
		value = (long) image[i][j];
	    }

	    byte = (unsigned char) (value % 256);
	    fwrite(&byte, sizeof(unsigned char), 1, out_file);

	    byte = (unsigned char) (value / 256);
	    fwrite(&byte, sizeof(unsigned char), 1, out_file);

	}
    }

    fclose(out_file);

    return;
} /* image_write_pgm_16bpp */

/*--------------------------------------------------------------------------*/

/** \brief Reads a 2D grey value image from a file in PDM format.
 */
void image_read_pdm(
    const char *file_name,   /**< name of file to read                    */
    double ***image, /**< the pointer to the image is stored here         */
    long   *nx,      /**< image size in x direction as given in the input 
                          file                                            */
    long   *ny,      /**< image size in y direction as given in the input 
                          file                                            */
    long   bx,       /**< boundary size in x direction                    */
    long   by,       /**< boundary size in y direction                    */
    struct comment **comments) /**< the comment field of the input file   */
{
    FILE* in_file;   /* input file handle   */
    char  row[256];  /* for reading data    */
    long  i, j;      /* loop variables      */
    long  value;     /* for data conversion */
    enum  pdm_storage_type st;  /* storage type F, I, B */

    /* open input file */
    in_file = fopen(file_name, "r");
    if(in_file == 0)
    {
	printf("image_read_pdm: Cannot open input file %s.\n", 
	       file_name);
	exit(1);
    }

    /* check file type with magic number */
    fgets(row, 255, in_file);
    if(strncmp(row, "P9", 2))
    {
	printf("image_read_pdm: File format is not PDM.\n");
	exit(1);
    }
    
    /* read comments */
    comment_alloc(comments);
    fgets(row, 255, in_file);
    while(row[0] == '#')
    {
	comment_add(*comments, row);
	fgets(row, 255, in_file);
        /*	printf("comment_read: %s", row);   */
    }

    if(strncmp(row, "SS", 2))
    {
	printf("image_read_pdm: Ordering is not correct.\n");
	exit(1);
    }

    /* read width (nx) and height (ny) */
    fgets(row, 255, in_file);
    sscanf(row, "%ld %ld", nx, ny);
    /* printf("nx: %ld, ny: %ld\n", *nx, *ny); */

    /* allocate memory for the image */
    image_alloc(image, *nx, *ny, bx, by);

    /* read storage type */
    fgets(row, 255, in_file);
    if(!strncmp(row, "F", 1))
    {
	st = F;
    } 
    else if(!strncmp(row, "I", 1))
    {
	st = I;
    } 
    else if(!strncmp(row, "I", 1))
    {
	st = B;
    }
    else 
    {
	printf("image_read_pdm: storage type %s is not supported.\n",
	       row);
	exit(1);
    }

    /* read image data */
    switch(st)
    {
	case F:
	    for(j=by; j<*ny+by; j++)
	    {
		for(i=bx; i<*nx+bx; i++)
		{
		    fscanf(in_file, "%lf", &((*image)[i][j]) );
		}
	    }
	    break;

	case I:
	    for(j=by; j<*ny+by; j++)
	    {
		for(i=bx; i<*nx+bx; i++)
		{
		    fscanf(in_file, "%ld", &value );
		    (*image)[i][j] = (double)value;
		}
	    }
	    break;
	    
	case B:
	    for(j=by; j<*ny+by; j++)
	    {
		for(i=bx; i<*nx+bx; i++)
		{
		    (*image)[i][j] = (double) getc(in_file);
		}
	    }
	    break;

	default:
	    printf("image_read_pdm: storage type %ld not supported\n",
		   (long) st);
	    exit(1);
    }

    fclose(in_file);
    printf("image_read_pdm: input file %s read.\n", file_name);

    return;
} /* image_read_pdm */

/*--------------------------------------------------------------------------*/

/** \brief Writes a 2D grey value image to a file in PDM format.
 */
void image_write_pdm(
    const char *file_name,      /**< name of file to write        */
    double **image,             /**< image data                   */
    long   nx,                  /**< image size in x direction    */
    long   ny,                  /**< image size in y direction    */
    long   bx,                  /**< boundary size in x direction */
    long   by,                  /**< boundary size in y direction */
    enum   pdm_storage_type st, /**< storage type for PDM file    */
    struct comment *comments)   /**< comment field                */
{
    FILE* out_file;
    long  i, j;
    unsigned char byte;

    /* open output file */
    out_file = fopen(file_name, "w");
    if(out_file == 0)
    {
	printf("image_write_pgm: Could not open output file %s.\n", 
	       file_name);
	exit(1);
    }

    /* ---- PGM header ---- */
    fprintf(out_file, "P9\n");

    /* print out comments */
    while(comments != NULL)
    {
	if(comments->line != NULL) 
	{
	    fprintf(out_file, "%s", comments->line);
	}
	comments = comments->next;
    }

    /* specify data format: */
    /* two space dimensions and only one grey value per pixel */
    fprintf(out_file, "SS\n");

    /* write image size in x and y direction */
    fprintf(out_file, "%ld %ld\n", nx, ny);

    /* write storage type */
    switch(st)
    {
	case F:
	    fprintf(out_file, "F\n");
	    break;
	
	case I:
	    fprintf(out_file, "I\n");
	    break;

	case B:
	    fprintf(out_file, "B\n");
	    break;

	default:
	    printf("image_write_pdm: storage type %ld not supported\n",
		   (long) st);
	    exit(1);
    }

    /* write image data */
    switch(st)
    {
	case F:
	    for(j=by; j<ny+by; j++)
	    {
		for(i=bx; i<nx+bx; i++)
		{
		    fprintf(out_file, "%e\n", image[i][j]);
		}
	    }
	    break;

	case I:
	    for(j=by; j<ny+by; j++)
	    {
		for(i=bx; i<nx+bx; i++)
		{
		    fprintf(out_file, "%ld\n", (long)image[i][j]);
		}
	    }
	    break;

	case B:
	    for(j=by; j<ny+by; j++)
	    {
		for(i=bx; i<nx+bx; i++)
		{
		    if(image[i][j] < 0) 
		    {
			byte = (unsigned char) 0.0;
		    } 
		    else if(image[i][j] > 255.0)
		    {
			byte = (unsigned char) 255.0;
		    } 
		    else
		    {
			byte = (unsigned char) image[i][j];
		    }
		    
		    fwrite(&byte, sizeof(unsigned char), 1, out_file);
		}
	    }
	    break;

	default:
	    printf("image_write_pdm: storage type %ld not supported\n",
		   (long) st);
	    exit(1);
    }

    fclose(out_file);

    return;
} /* image_write_pdm */

/*--------------------------------------------------------------------------*/

#endif /* __IMIO_H__ */

