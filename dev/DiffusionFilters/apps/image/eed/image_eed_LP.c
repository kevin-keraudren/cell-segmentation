/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                        EDGE-ENHANCING DIFFUSION                          */
/*                                                                          */
/* J. Weickert,                                                             */
/* Theoretical foundations of anisotropic diffusion in image processing,    */
/* Computing, Suppl. 11, 221-236, 1996.                                     */
/*                                                                          */
/* Explicit time discretisation.                                            */
/*                                                                          */
/* Maintainer: Stephan Didas <didas@mia.uni-saarland.de>                    */
/*             Luis Pizarro <luis.pizarro@imperial.ac.uk>                   */
/*                                                                          */
/*--------------------------------------------------------------------------*/
/* Ad-hoc stuff for LP:                                                     */
/*                                                                          */
/* + all parameters passed in the command line                              */
/* + new 'noOutImages' parameter: how many output images will be written    */
/* + 'out_file' is used as prefix for the output images 'out_file_no*.pgm'  */
/*                                                                          */
/* + removed 'image_ced_ex', calling directly 'image_ced_ex_step'           */
/*--------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "base/memory.h"
#include "image/imana.h"
#include "image/imio.h"
#include "image/immem.h"
#include "image/imconv.h"
#include "image/imst.h"
#include "matrix/matrix.h"
#include "matrix/mat_sym2.h"
#include "tensor2d/ten2dmem.h"

/*--------------------------------------------------------------------------*/

/** \brief Determines how the weight \f$theta\f$ is chosen in the stencil 
 for anisotropic nonlinear diffusion filtering. 
 */
enum aniso_discr {
	MINUS_ONE,  /**< \f$\theta = -1$ in all pixels */
	ZERO,       /**< \f$\theta =  0$ in all pixels */
	ONE,        /**< \f$\theta = +1$ in all pixels */
	ADAPTIVE    /**< \f$\theta = \sign (d_{1,0} (i, j))\f$ */
};

/*--------------------------------------------------------------------------*/

/** \brief Performs one explicit step of edge-enhancing diffusion 
 for grey value images.
 */
void image_eed_ex_step(
											 double **f,     /**< input: original image, unchanged       */  
											 long   nx,      /**< image size in x direction              */
											 long   ny,      /**< image size in y direction              */
											 long   bx,      /**< boundary size in x direction           */
											 long   by,      /**< boundary size in y direction           */
											 double hx,      /**< spatial step size in x direction       */
											 double hy,      /**< spatial step size in y direction       */
											 double sigma,   /**< noise scale for structure tensor       */
											 enum fd_scheme scheme, /**< derivative approximation scheme */
											 double prec,    /**< cutoff precision for convolutions      */
											 double cparam,  /**< contrast parameter lambda              */
											 enum aniso_discr discr, /**< how to set theta               */
											 double tau,     /**< time step size                         */
											 double ***dt,   /**< diffusion tensor                       */
											 double **theta, /**< discretisation parameter               */
											 double ***stencil,   /**< entries of the evolution matrix   */
											 double **u)     /**< output image                           */
{
	long   i, j;      /* loop variables */
	
	double **d;       /* diffusion tensor in matrix notation */
	double *lambda;   /* eigenvalues of diffusion tensor     */
	double **v;       /* eigenvectors of diffusion tensor    */
	
	double delta = 10e-14;   /* machine precision for PA_trans */
	double eps   = 10e-06;   /* error bound for PA_trans       */
	
	/* factors for anisotropic diffuion stencil */
	double hxx_den = 1.0 / (2.0*hx*hx);
	double hyy_den = 1.0 / (2.0*hy*hy);
	double hxy_den = 1.0 / (4.0*hx*hy);
	
	double oocontr_param = 1.0 / (cparam*cparam);
	
	if( bx < 1 || by < 1 )
	{
		printf("image_eed_LP: boundary size too small.\n");
		exit(1);
	}
	
	matrix_alloc(&d, 4, 4);
	vector_alloc(&lambda, 4);
	matrix_alloc(&v, 4, 4);
	
	/* calculate the structure tensors */
	image_structure_tensor(f, nx, ny, bx, by, hx, hy, sigma, scheme,
												 0.0, prec, dt);
	
	/* calculate the diffusion tensors */
	for(i=bx; i<nx+bx; i++)
	{
		for(j=bx; j<ny+by; j++)
		{
	    /* calculate principal axis transformation */
	    mat_sym2_to_matrix(dt[i][j], d);
	    PA_trans(d, 2, delta, eps, lambda, v);
			
	    /* transform eigenvalues of the tensor */
	    lambda[1] = 1.0 / (1.0 + (lambda[1]*lambda[1]) * oocontr_param); 
			lambda[2] = 1.0;
			
	    /* backtransform of the tensor and writing back */
	    PA_backtrans(lambda, v, 2, d);
	    mat_sym2_from_matrix(d, dt[i][j]);
		}
	}
	ten2d_mirror_bd(dt, nx, ny, bx, by, 3);
	
	/* set discretisation parameter theta */
	switch(discr) {
		case MINUS_ONE:
	    for(i=bx; i<nx+bx; i++)
	    {
				for(j=by; j<ny+by; j++)
				{
					theta[i][j] = -1.0;
				}
	    }
	    break;
	    
		case ZERO:
	    for(i=bx; i<nx+bx; i++)
	    {
				for(j=by; j<ny+by; j++)
				{
					theta[i][j] = 0.0;
				}
	    }
	    break;
			
		case ONE:
	    for(i=bx; i<nx+bx; i++)
	    {
				for(j=by; j<ny+by; j++)
				{
					theta[i][j] = 1.0;
				}
	    }
	    break;
			
		case ADAPTIVE:
	    for(i=bx; i<nx+bx; i++)
	    {
				for(j=by; j<ny+by; j++)
				{
					if( dt[i][j][2] > 0.0 ) 
					{
						theta[i][j] = 1.0;
					} 
					else if( dt[i][j][2] < 0.0 )
					{
						theta[i][j] = -1.0;
					} 
					else 
					{
						theta[i][j] = 0.0;
					}
				}
	    }
	    break;
	}
	
	/* calculate the explicit step */
	image_mirror_bd(f, nx, ny, bx, by);
	
	for(i=bx; i<nx+bx; i++)
	{
		for(j=by; j<ny+by; j++)
		{
	    stencil[i][j][0] = hxy_den * (i>bx) * (j>by) 
			* ( + (1+theta[i-1][j]) * dt[i-1][j][2]
				 + (1+theta[i][j-1]) * dt[i][j-1][2] );
			
	    stencil[i][j][1] = hxy_den * ( 
																		+ (i>bx) * (j>by) * ( 
																												 + (1-theta[i  ][j]) * dt[i  ][j][2]
																												 - (1+theta[i-1][j]) * dt[i-1][j][2] )
																		+ (i>bx) * (j<ny+by-1) * (
																															- (1+theta[i  ][j]) * dt[i  ][j][2]
																															+ (1-theta[i-1][j]) * dt[i-1][j][2] ) )
			+ hxx_den * (i>bx) * (dt[i][j][0] + dt[i-1][j][0]);
			
	    stencil[i][j][2] = hxy_den * (i>bx) * (j<ny+by-1)
			* ( - (1-theta[i-1][j]) * dt[i-1][j][2] 
				 - (1-theta[i][j+1]) * dt[i][j+1][2] );
			
	    stencil[i][j][3] = hxy_den * (
																		+ (i<nx+bx-1) * (j>by) * (
																															+ (1-theta[i][j-1]) * dt[i][j-1][2]
																															- (1+theta[i][j  ]) * dt[i][j  ][2] )
																		+ (i>bx) * (j>by) * (
																												 + (1-theta[i][j  ]) * dt[i][j  ][2] 
																												 - (1+theta[i][j-1]) * dt[i][j-1][2] ) )
			+ (j>by) * hyy_den * ( dt[i][j][1] + dt[i][j-1][1] );
		}
	}
	
	/* set boundaries to zero */
	for(i=bx; i<nx+bx; i++)
	{
		for(j=0; j<by; j++)
		{
	    stencil[i][j][0] = 0.0;
	    stencil[i][j][1] = 0.0;
	    stencil[i][j][2] = 0.0;
	    stencil[i][j][3] = 0.0;
		}
		
		for(j=ny+by; j<ny+2*by; j++)
		{
	    stencil[i][j][0] = 0.0;
	    stencil[i][j][1] = 0.0;
	    stencil[i][j][2] = 0.0;
	    stencil[i][j][3] = 0.0;
		}
	}
	
	/* set boundaries to zero */
	for(j=0; j<ny+2*by; j++)
	{
		for(i=0; i<bx; i++)
		{
	    stencil[i][j][0] = 0.0;
	    stencil[i][j][1] = 0.0;
	    stencil[i][j][2] = 0.0;
	    stencil[i][j][3] = 0.0;
		}
		
		for(i=nx+bx; i<nx+2*bx; i++)
		{
	    stencil[i][j][0] = 0.0;
	    stencil[i][j][1] = 0.0;
	    stencil[i][j][2] = 0.0;
	    stencil[i][j][3] = 0.0;
		}
	}
	
	/* derive missing matrix entries by symmetry */
	for(i=bx; i<nx+bx; i++)
	{
		for(j=by; j<ny+by; j++)
		{
	    stencil[i][j][5] = stencil[i  ][j+1][3];
	    stencil[i][j][6] = stencil[i+1][j-1][2];
	    stencil[i][j][7] = stencil[i+1][j  ][1];
	    stencil[i][j][8] = stencil[i+1][j+1][0];
			
	    /* central element */
	    stencil[i][j][4] = 
			- stencil[i][j][0] - stencil[i][j][1]
			- stencil[i][j][2] - stencil[i][j][3]
			- stencil[i][j][5] - stencil[i][j][6]
			- stencil[i][j][7] - stencil[i][j][8] ;
		}
	}    
	
	for(i=bx; i<nx+bx; i++)
	{
		for(j=by; j<ny+by; j++)
		{
	    /* explicit step */
	    u[i][j] = f[i][j] + tau *
			(   stencil[i][j][0] * f[i-1][j-1]
			 + stencil[i][j][1] * f[i-1][j  ]
			 + stencil[i][j][2] * f[i-1][j+1]
			 + stencil[i][j][3] * f[i  ][j-1]
			 + stencil[i][j][4] * f[i  ][j  ]
			 + stencil[i][j][5] * f[i  ][j+1]
			 + stencil[i][j][6] * f[i+1][j-1]
			 + stencil[i][j][7] * f[i+1][j  ]
			 + stencil[i][j][8] * f[i+1][j+1] );
		}
	}
	
	matrix_disalloc(d);
	vector_disalloc(lambda);
	matrix_disalloc(v);
	
	return;
} /* image_eed_ex_step */

/*--------------------------------------------------------------------------*/

int main(int argc, char* argv[])
{
	char   in_file[256];   /* input file name */
	char   out_file[256];  /* output file name */
	char   prefix_out_file[256];  /* prefix file name */
  
	double **f;            /* initial image data */
	double **u;            /* iteration result */
	
	long   nx, ny;         /* image size in x, y direction */
	long   bx = 1;         /* boundary size in x direction */
	long   by = 1;         /* boundary size in y direction */
	
	double sigma;          /* noise scale for structure tensor       */
	int    s1, s2;         /* for user inputs                        */
	enum fd_scheme scheme; /* derivative approximation scheme        */
	double prec;           /* cutoff precision for convolutions      */
	
	double cparam;          /* contrast parameter for eed */
	enum aniso_discr discr; /* discretisation parameter theta */
	
	double tau;             /* time step size */
	long   it, it_step;     /* number of iterations */
	long   k;
	
	struct comment* comments;     /* comments in the input file  */
	struct comment* new_comments; /* new comment for output file */
	
	double ***dt;   /* diffusion tensor */
	double ***stencil;   /* entries of the evolution matrix */
	double **theta; /* discretisation parameter */
	
	if (argc < 6) {
		printf("\nError: you must provide 5 parameters in the command line.\n");
		printf("Usage:\n");
		printf("image_eed_LP [input_image.PGM] [cparam] [iters_total] [iters_step] [prefix_out_image (without .PGM)]\n\n");
		exit(0);
	}

	/* read parameters */
	strcpy(in_file,argv[1]);
	sigma = 0.5f;
	s1 = 0; /* scheme = 0 */
	cparam = atof(argv[2]);
	prec = 3.0f;
	s2 = 0; /* theta = 0 */
	tau = 0.2f;
	it = atol(argv[3]);
	it_step = atol(argv[4]);
	strcpy(prefix_out_file,argv[5]);
	
	/* print out parameters */
	printf("parameters:\n");	
  printf("input_image = %s\n", in_file);
  printf("sigma = %f\n", sigma);
  printf("scheme = %d\n", s1);
  printf("cparam = %f\n", cparam);
  printf("prec = %f\n", prec);
  printf("theta = %d\n", s2);
  printf("tau = %f\n", tau);
  printf("iters_total = %ld\n", it);
  printf("iters_step = %ld\n", it_step);
	printf("prefix_out_image = %s\n", prefix_out_file);
	
	switch(s1) 
	{
		case 0:
	    scheme = CENTRAL;
	    break;
			
		case 1:
	    scheme = SOBEL;
	    break;
			
		case 2:
	    scheme = KUMAR;
	    break;
			
		default:
	    fprintf(stderr, "image_eed_LP: scheme must be 0, 1 or 2.");
	    exit(1);
	}
	
	switch(s2) 
	{
		case 0:
	    discr = ZERO;
	    break;
			
		case 1:
	    discr = ONE;
	    break;
			
		case 2:
	    discr = MINUS_ONE;
	    break;
			
		case 3:
	    discr = ADAPTIVE;
	    break;
			
		default:
	    fprintf(stderr, "image_eed_LP: theta must be 0, 1, 2 or 3.");
	    exit(1);
	}
	
	/* read image name */
	image_read_pgm(in_file, &f, &nx, &ny, bx, by, &comments);
	image_alloc(&u, nx, ny, bx, by);
	
	/* prepare comments for output file */
	comment_alloc(&new_comments);
	comment_add_sep(new_comments);
	comment_add(new_comments, "# %s\n", out_file);
	comment_add(new_comments, "# initial image:  %s\n", in_file);
	comment_add(new_comments, "# Edge enhancing diffusion filtering\n");
	comment_add(new_comments, "# structure tensor parameters:\n");
	comment_add(new_comments, "#   noise scale sigma:     %f\n", sigma);
	switch(scheme) 
	{
		case CENTRAL:
	    comment_add(new_comments, "#   central differences\n");
	    break;
			
		case 1:
	    comment_add(new_comments, "#   Sobel operators\n");
	    break;
			
		case 2:
	    comment_add(new_comments, "#   Kumar operators\n");
	    break;
			
		default:
	    fprintf(stderr, "image_eed_LP: finite difference scheme %d "
							"not known\n", s1);
	    exit(1);
	}
	comment_add(new_comments, "#   cutoff precision:      %f\n", prec);
	comment_add(new_comments, "# EED contrast parameter:  %f\n", cparam);
	comment_add(new_comments, "# discretisation parameter theta:\n");
	switch(discr) 
	{
		case ZERO:
	    comment_add(new_comments, "#   theta = 0\n");
	    break;
			
		case ONE:
	    comment_add(new_comments, "#   theta = +1\n");
	    break;
			
		case MINUS_ONE:
	    comment_add(new_comments, "#   theta = -1\n");
	    break;
			
		case ADAPTIVE:
	    comment_add(new_comments, "#   theta = sign( d_{1,0}) \n");
	    break;
			
		default:
	    fprintf(stderr, "image_eed_LP: discretisation parameter %d "
							"not known\n", discr);
	    exit(1);
	}
	comment_add(new_comments, "# time step size:          %f\n", tau);
	comment_add(new_comments, "# iteration number:        %ld\n", it);
	comment_add(new_comments, "# iteration step:          %ld\n", it_step);
	comment_add(new_comments, "# --> stopping time:       %f\n",
							(double)it * tau);
	comment_add(new_comments, "# Created with %s\n", argv[0]);
	comment_add_time(new_comments);
	comment_add_sep(new_comments);
	comment_append(new_comments, comments);
	
	if( bx < 1 || by < 1 )
	{
		printf("image_eed_LP: boundary size too small.\n");
		exit(1);
	}
	
	ten2d_alloc(&dt, nx, ny, bx, by, 3);
	array3_alloc(&stencil, nx+2*bx, ny+2*by, 9);
	matrix_alloc(&theta, nx+2*bx, ny+2*by);
	
	for(k=1; k<=it; k++)
	{
		
		image_eed_ex_step(f, nx, ny, bx, by, 1.0, 1.0, sigma, scheme, prec, 
											cparam, discr, tau, dt, theta, stencil, u);
				
		/* write multiple output files according to it_steps */
		if (k % it_step == 0) {
			printf("iteration number %5ld:\n", k);
	    image_analyse_print(u, nx, ny, bx, by, "current image");
			
			sprintf(out_file, "%s_%04ld.pgm",prefix_out_file,k);
			image_write_pgm(out_file, u, nx, ny, bx, by, new_comments);
			printf("output image written: %s\n", out_file);
		}
		
		/* copy result back to f */
		image_copy(u, nx, ny, bx, by, f);
		
	}
		
	/* ---- disallocate storage ---- */
	image_disalloc(f);
	image_disalloc(u);
	comment_disalloc(new_comments);
	
	ten2d_disalloc(dt);
	array3_disalloc(stencil);
	matrix_disalloc(theta);

	return(0);
} /* main */

/*--------------------------------------------------------------------------*/
