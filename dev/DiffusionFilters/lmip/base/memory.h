/*--------------------------------------------------------------------------*/
/** \file memory.h

   \brief Memory allocation methods for 1D double vectors.

   Maintainer: Stephan Didas <didas@mia.uni-saarland.de>                    */
/*--------------------------------------------------------------------------*/

#ifndef __MEMORY_H__
#define __MEMORY_H__

#include <stdlib.h>
#include <stdio.h>

/*--------------------------------------------------------------------------*/

/** \brief Allocates storage for a vector of size nx.

    Exits the program with value 1 if there is not enough storage available.
 */
void vector_alloc(double **vector,      /**< the pointer to the vector is 
                                             stored here                   */
                  long   nx)            /**< size of the vector            */
{
    *vector = (double *) malloc (nx * sizeof(double));
    
    if (*vector == NULL)
    {
	printf("vector_alloc: not enough storage available\n");
	exit(1);
    }
    
    return;
} /* vector_alloc */

/*--------------------------------------------------------------------------*/

/** Disallocates storage of a vector. 
 */
void vector_disalloc(double *vector)    /**< pointer to the vector  */
{
    free(vector);

    return;
} /* vector_disalloc */

/*--------------------------------------------------------------------------*/

/** \brief Allocates storage for a matrix of size nx * ny.

    The storage for the entries of the matrix is allocated as one linear 
    memory block. 
 */
void matrix_alloc(double ***matrix, /**< pointer to the matrix */
		  long   nx,        /**< size in x direction   */
		  long   ny)        /**< size in y direction   */
{
    long   i;     /* loop variable */
    double *tmp;  /* time saver */
    
    tmp = (double*) malloc (nx * ny * sizeof(double));
    if (tmp == NULL)
    {
	printf("matrix_alloc: not enough storage available\n");
	exit(1);
    }
    
    *matrix = (double**) malloc (nx * sizeof(double*));
    if (*matrix == NULL)
    {
	printf("alloc_matrix: not enough storage available\n");
	exit(1);
    }
    
    for (i=0; i<nx; i++)
    {
	(*matrix)[i] = &(tmp[ny*i]);
    }

    return;
} /* matrix_alloc */

/*--------------------------------------------------------------------------*/

/** \brief Disallocates storage of a matrix of size nx * ny.
 */
void matrix_disalloc(double **matrix) /**< the matrix */
{
    free(matrix[0]);
    free(matrix);

    return;
} /* matrix_disalloc */

/*--------------------------------------------------------------------------*/

/** \brief Allocates storage for a three-dimensional array of size 
    n1*n2*n3.

    Exits the program with value 1 if there is not enough storage 
    available.
 */
void array3_alloc(double ****array,   /**< stores the pointer to the array  */
 		  long   n1,          /**< array size in direction 1        */
		  long   n2,          /**< array size in direction 2        */
		  long   n3)          /**< array size in direction 3        */
{
    long i, j;

    double *tmp01;   /* the data storage     */
    double **tmp02;  /* first pointer array  */

    tmp01 = (double*) malloc( n1 * n2 * n3 * sizeof(double) );
    if(tmp01 == NULL) 
    {
	printf("array3_alloc: not enough storage available for tmp01.\n");
	exit(1);
    }

    tmp02 = (double**) malloc( n1 * n2 * sizeof(double*) );
    if(tmp02 == NULL)
    {
	printf("array3_alloc: not enough storage available for tmp02.\n");
	exit(1);
    }

    *array = (double***) malloc( n1 * sizeof(double**) );
    if(*array == NULL)
    {
	printf("array3_alloc: not enough storage available for *array.\n");
	exit(1);
    }

    for(i=0; i<n1; i++)
    {
        for(j=0; j<n2; j++)
        {
            tmp02[i*n2 + j] = &( tmp01[ (i*n2 + j)*n3 ] );
        }

        (*array)[i] = &( tmp02[i*n2] );
    }
    
    return;
} /* array3_alloc */

/*--------------------------------------------------------------------------*/

/** \brief Disallocates storage of a three-dimensional array.
 */
void array3_disalloc(double ***array)   /**< the array */
{
    free(array[0][0]);
    free(array[0]);
    free(array);

    return;
} /* array3_disalloc */

/*--------------------------------------------------------------------------*/

/** \brief Allocates storage for a four-dimensional array of size 
    n1*n2*n3*n4.

    Exits the program with value 1 if there is not enough storage 
    available.
 */
void array4_alloc(double *****array,  /**< stores the pointer to the array  */
 		  long   n1,          /**< array size in direction 1        */
		  long   n2,          /**< array size in direction 2        */
		  long   n3,          /**< array size in direction 3        */
		  long   n4)          /**< array size in direction 4        */
{
    long i, j, k;

    double *tmp01;   /* the data storage     */
    double **tmp02;  /* first pointer array  */
    double ***tmp03; /* second pointer array */

    tmp01 = (double*) malloc( n1 * n2 * n3 * n4 * sizeof(double) );
    if(tmp01 == NULL)
    {
	printf("array4_alloc: not enough storage available for tmp01.\n");
	exit(1);
    }

    tmp02 = (double**) malloc( n1 * n2 * n3 * sizeof(double*) );
    if(tmp02 == NULL)
    {
	printf("array4_alloc: not enough storage available for tmp02.\n");
	exit(1);
    }

    tmp03 = (double***) malloc( n1 * n2 * sizeof(double**) );
    if(tmp03 == NULL)
    {
	printf("array4_alloc: not enough storage available for tmp03.\n");
	exit(1);
    }

    *array = (double****) malloc( n1 * sizeof(double***) );
    if(*array == NULL)
    {
	printf("array4_alloc: not enough storage available for *array.\n");
	exit(1);
    }

    for(i=0; i<n1; i++)
    {
        for(j=0; j<n2; j++)
        {
            for(k=0; k<n3; k++)
            {
                tmp02[(i*n2 + j)*n3 + k] = 
                    &( tmp01[( (i*n2 + j)*n3 + k )*n4] );
            }

            tmp03[i*n2 + j] = 
                  &( tmp02[ (i*n2 + j)*n3 ] );
        }

        (*array)[i] = &( tmp03[i*n2] );
    }
    
    return;
} /* array4_alloc */

/*--------------------------------------------------------------------------*/

/** \brief Disallocates storage of a four-dimensional array.
 */
void array4_disalloc(double ****array)   /**< the array */
{
    free(array[0][0][0]);
    free(array[0][0]);
    free(array[0]);
    free(array);

    return;
} /* array4_disalloc */

/*--------------------------------------------------------------------------*/

/** \brief Allocates storage for a five-dimensional array of size 
    n1*n2*n3*n4*n5.

    Exits the program with value 1 if there is not enough storage 
    available.
 */
void array5_alloc(double ******array, /**< stores the pointer to the array  */
 		  long   n1,          /**< array size in direction 1        */
		  long   n2,          /**< array size in direction 2        */
		  long   n3,          /**< array size in direction 3        */
		  long   n4,          /**< array size in direction 4        */
		  long   n5)          /**< array size in direction 5        */
{
    long i, j, k, l;

    double *tmp01;    /* the data storage     */
    double **tmp02;   /* first pointer array  */
    double ***tmp03;  /* second pointer array */
    double ****tmp04; /* third pointer array  */

    tmp01 = (double*) malloc( n1 * n2 * n3 * n4 * n5 * sizeof(double) );
    if(tmp01 == NULL)
    {
	printf("array5_alloc: not enough storage available for tmp01.\n");
	exit(1);
    }

    tmp02 = (double**) malloc( n1 * n2 * n3 * n4 * sizeof(double*) );
    if(tmp02 == NULL)
    {
	printf("array5_alloc: not enough storage available for tmp02.\n");
	exit(1);
    }

    tmp03 = (double***) malloc( n1 * n2 * n3 * sizeof(double**) );
    if(tmp03 == NULL)
    {
	printf("array5_alloc: not enough storage available for tmp03.\n");
	exit(1);
    }

    tmp04 = (double****) malloc( n1 * n2 * sizeof(double***) );
    if(tmp04 == NULL)
    {
	printf("array5_alloc: not enough storage available for tmp04.\n");
	exit(1);
    }

    *array = (double*****) malloc( n1 * sizeof(double****) );
    if(array == NULL)
    {
	printf("array5_alloc: not enough storage available for array.\n");
	exit(1);
    }

    for(i=0; i<n1; i++)
    {
        for(j=0; j<n2; j++)
        {
            for(k=0; k<n3; k++)
            {
		for(l=0; l<n4; l++)
		{
		    tmp02[((i*n2 + j)*n3 + k)*n4 + l] =
			&( tmp01[( ( ( i*n2 + j )*n3 + k)*n4 + l)*n5] );
		}

                tmp03[(i*n2 + j)*n3 + k] = 
                    &( tmp02[( (i*n2 + j)*n3 + k )*n4] );
            }

            tmp04[i*n2 + j] = 
                  &( tmp03[ (i*n2 + j)*n3 ] );
        }

        (*array)[i] = &( tmp04[i*n2] );
    }
    
    return;
} /* array5_alloc */

/*--------------------------------------------------------------------------*/

/** \brief Disallocates storage of a five-dimensional array.
 */
void array5_disalloc(double *****array)   /**< the array */
{
    free(array[0][0][0][0]);
    free(array[0][0][0]);
    free(array[0][0]);
    free(array[0]);
    free(array);

    return;
} /* array5_disalloc */

/*--------------------------------------------------------------------------*/

/** \brief Allocates storage for a matrix of size nx * ny of chars.

    The storage for the entries of the matrix is allocated as one linear 
    memory block. 
 */
void char_matrix_alloc(char ***matrix, /**< pointer to the matrix */
		       long   nx,        /**< size in x direction   */
		       long   ny)        /**< size in y direction   */
{
    long   i;     /* loop variable */
    char *tmp;    /* time saver */
    
    tmp = (char*) malloc (nx * ny * sizeof(char));
    if (tmp == NULL)
    {
	printf("matrix_alloc: not enough storage available\n");
	exit(1);
    }
    
    *matrix = (char**) malloc (nx * sizeof(char*));
    if (*matrix == NULL)
    {
	printf("alloc_matrix: not enough storage available\n");
	exit(1);
    }
    
    for (i=0; i<nx; i++)
    {
	(*matrix)[i] = &(tmp[ny*i]);
    }

    return;
} /* matrix_alloc */

/*--------------------------------------------------------------------------*/

/** \brief Disallocates storage of a matrix of size nx * ny.
 */
void char_matrix_disalloc(char **matrix) /**< the matrix */
{
    free(matrix[0]);
    free(matrix);

    return;
} /* char_matrix_disalloc */

/*--------------------------------------------------------------------------*/

/** \brief Allocates storage for a three-dimensional array of chars of 
    size n1*n2*n3.

    Exits the program with value 1 if there is not enough storage 
    available.
 */
void char_array3_alloc(
    char ****array,   /**< stores the pointer to the array  */
    long n1,          /**< array size in direction 1        */
    long n2,          /**< array size in direction 2        */
    long n3)          /**< array size in direction 3        */
{
    long i, j;

    char *tmp01;   /* the data storage     */
    char **tmp02;  /* first pointer array  */

    tmp01 = (char*) malloc( n1 * n2 * n3 * sizeof(char) );
    if(tmp01 == NULL) 
    {
	printf("array3_alloc: not enough storage available for tmp01.\n");
	exit(1);
    }

    tmp02 = (char**) malloc( n1 * n2 * sizeof(char*) );
    if(tmp02 == NULL)
    {
	printf("array3_alloc: not enough storage available for tmp02.\n");
	exit(1);
    }

    *array = (char***) malloc( n1 * sizeof(char**) );
    if(*array == NULL)
    {
	printf("array3_alloc: not enough storage available for *array.\n");
	exit(1);
    }

    for(i=0; i<n1; i++)
    {
        for(j=0; j<n2; j++)
        {
            tmp02[i*n2 + j] = &( tmp01[ (i*n2 + j)*n3 ] );
        }

        (*array)[i] = &( tmp02[i*n2] );
    }
    
    return;
} /* char_array3_alloc */

/*--------------------------------------------------------------------------*/

/** \brief Disallocates storage of a three-dimensional array of chars.
 */
void char_array3_disalloc(char ***array)   /**< the array */
{
    free(array[0][0]);
    free(array[0]);
    free(array);

    return;
} /* char_array3_disalloc */

/*--------------------------------------------------------------------------*/

#endif /* __MEMORY_H__ */

