/*--------------------------------------------------------------------------*/
/** \file imana.h 

   \brief Analysis methods for 2D double images. 

   Maintainer: Stephan Didas <didas@mia.uni-saarland.de>                    */
/*--------------------------------------------------------------------------*/

#ifndef __IMANA_H__
#define __IMANA_H__

#include <math.h>

/*--------------------------------------------------------------------------*/

/** \brief Calculates the l1 norm of an image u.
 */
double image_l1norm(double **u, /**< the image, unchanged */
		    long   nx,  /**< image size in x direction    */
		    long   ny,  /**< image size in y direction    */
		    long   bx,  /**< boundary size in x direction */
		    long   by)  /**< boundary size in y direction */
{
    long i, j;
    double norm = 0.0;

    for(i=bx; i<nx+bx; i++)
    {
	for(j=by; j<ny+by; j++)
	{
	    norm = norm + fabs(u[i][j]);
	}
    }

    return norm;
} /* image_l1norm */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the l2 norm of an image u.
 */
double image_l2norm(double **u, /**< the image, unchanged         */
		    long   nx,  /**< image size in x direction    */
		    long   ny,  /**< image size in y direction    */
		    long   bx,  /**< boundary size in x direction */
		    long   by)  /**< boundary size in y direction */
{
    long i, j;
    double norm = 0.0;

    for(i=bx; i<nx+bx; i++)
    {
	for(j=by; j<ny+by; j++)
	{
	    norm = norm + u[i][j]*u[i][j];
	}
    }

    return sqrt(norm);
} /* image_l2norm */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the l_infinity norm of an image u.
 */
double image_linfnorm(double **u, /**< the image, unchanged        */
		      long   nx,  /**< image size in x direction    */
		      long   ny,  /**< image size in y direction    */
		      long   bx,  /**< boundary size in x direction */
		      long   by)  /**< boundary size in y direction */
{
    long i, j;
    double help;
    double norm = fabs(u[bx][by]);

    for(i=bx; i<nx+bx; i++)
    {
	for(j=by; j<ny+by; j++)
	{
	    help = fabs(u[i][j]);
	    if(help > norm) norm = help;
	}
    }

    return norm;
} /* image_linfnorm */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the lp norm of an image u.
 */
double image_lpnorm(double **u,       /**< the image, unchanged         */
		    long   nx,        /**< image size in x direction    */
		    long   ny,        /**< image size in y direction    */
		    long   bx,        /**< boundary size in x direction */
		    long   by,        /**< boundary size in y direction */
		    unsigned long p)  /**< index of the lp space        */
{
    long i, j;
    double norm = 0.0;

    for(i=bx; i<nx+bx; i++)
    {
	for(j=by; j<ny+by; j++)
	{
	    norm = norm + pow(fabs(u[i][j]), p);
	}
    }

    return pow(norm, 1.0 / (double)p);
} /* image_lpnorm */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the l1 norm of the difference of two images f and u.
 */
double image_diff_l1norm(double **f,   /**< first image, unchanged       */
			 double **u,   /**< second image, unchanged      */
			 long   nx,    /**< image size in x direction    */
			 long   ny,    /**< image size in y direction    */
			 long   bx,    /**< boundary size in x direction */
			 long   by)    /**< boundary size in y direction */
{
    long i, j;
    double norm = 0.0;

    for(i=bx; i<nx+bx; i++)
    {
	for(j=by; j<ny+by; j++)
	{
	    norm = norm + fabs(u[i][j]-f[i][j]);
	}
    }

    return norm;
} /* image_diff_l1norm */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the l2 norm of the difference of two images f and u.
 */
double image_diff_l2norm(double **f,   /**< first image, unchanged       */
			 double **u,   /**< second image, unchanged      */
			 long   nx,    /**< image size in x direction    */
			 long   ny,    /**< image size in y direction    */
			 long   bx,    /**< boundary size in x direction */
			 long   by)    /**< boundary size in y direction */
{
    long i, j;
    double norm = 0.0;

    for(i=bx; i<nx+bx; i++)
    {
	for(j=by; j<ny+by; j++)
	{
	    norm = norm + (u[i][j]-f[i][j])*(u[i][j]-f[i][j]);
	}
    }

    return sqrt(norm);
} /* image_diff_l2norm */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the averaged mean square error of an image.
 */
double image_amse(double **f,   /**< first image, unchanged       */
		  double **u,   /**< second image, unchanged      */
		  long   nx,    /**< image size in x direction    */
		  long   ny,    /**< image size in y direction    */
		  long   bx,    /**< boundary size in x direction */
		  long   by)    /**< boundary size in y direction */
{
    long i, j;
    double norm = 0.0;

    for(i=bx; i<nx+bx; i++)
    {
	for(j=by; j<ny+by; j++)
	{
	    norm = norm + (u[i][j]-f[i][j])*(u[i][j]-f[i][j]);
	}
    }

    return norm / (double)(nx*ny);
} /* image_amse */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the l_infinity norm of the difference of two images 
           f and u.
 */
double image_diff_linfnorm(double **f,  /**< first image, unchanged       */
			   double **u,  /**< second image, unchanged      */
			   long   nx,   /**< image size in x direction    */
			   long   ny,   /**< image size in x direction    */
			   long   bx,   /**< boundary size in x direction */
			   long   by)   /**< boundary size in x direction */
{
    long i, j;
    double help;
    double norm = fabs(u[bx][by]-f[bx][by]);

    for(i=bx; i<nx+bx; i++)
    {
	for(j=by; j<ny+by; j++)
	{
	    help = fabs(u[i][j]-f[i][j]);
	    if(help > norm) norm = help;
	}
    }

    return norm;
} /* image_diff_linfnorm */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the lp norm of the difference of two images f and u.
 */
double image_diff_lpnorm(double **f,      /**< first image, unchanged  */
			 double **u,      /**< second image, unchagend */
			 long   nx,       /**< image size in x direction    */
			 long   ny,       /**< image size in y direction    */
			 long   bx,       /**< boundary size in x direction */
			 long   by,       /**< boundary size in y direction */
			 unsigned long p) /**< index of the lp space    */
{
    long i, j;
    double norm = 0.0;

    for(i=bx; i<nx+bx; i++)
    {
	for(j=by; j<ny+by; j++)
	{
	    norm = norm + pow(fabs(u[i][j]-f[i][j]), p);
	}
    }

    return pow(norm, 1.0 / (double)p);
} /* image_diff_lpnorm */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the Signal-to-Noise Ratio (SNR) for a degraded 
           version with respect to its original image.

    Given an original image \f$g\f$ and a degraded version \f$f\f$, the 
    Signal-to-Noise Ratio (SNR) is defined as
    \f[
      \mbox{SNR}(f, g) := 10 \log_{10} 
      \left( 
        \frac{\sum\limits_{i=1}^{nx} \sum\limits_{j=1}^{ny} 
	      (g_{ij} - \mu)^2}
	     {\sum\limits_{i=1}^{nx} \sum\limits_{j=1}^{ny} 
	      (f_{ij} - g_{ij})^2}
      \right) \enspace ,
    \f]
    where \f$\mu := \frac{1}{nx \cdot ny} 
    \sum\limits_{i=1}^{nx} \sum\limits_{j=1}^{ny} g_{ij}\f$ 
    denotes the mean value of the original image \f$g\f$.

    Source: Lecture Notes to 
    'Image Processing and Computer Vision', winter term 2005/2006,
    Prof. Weickert (Lecture 02, Slide 12-13)
 */
double image_snr(double **g,   /**< original image, unchanged  */
		 double **f,   /**< degraded image, unchagend */
		 long   nx,    /**< image size in x direction    */
		 long   ny,    /**< image size in y direction    */
		 long   bx,    /**< boundary size in x direction */
		 long   by)    /**< boundary size in y direction */
{
    long i, j;
    double mean = 0.0;
    double num  = 0.0;
    double den  = 0.0;

    /* determine the mean of the original image g */
    for(i=bx; i<nx+bx; i++)
    {
	for(j=by; j<ny+by; j++)
	{
	    mean = mean + g[i][j];
	}
    }
    mean = mean / (double)(nx * ny);

    for(i=bx; i<nx+bx; i++)
    {
	for(j=by; j<ny+by; j++)
	{
	    num = num + ( g[i][j] - mean ) * ( g[i][j] - mean );
	    den = den + ( f[i][j] - g[i][j] ) * ( f[i][j] - g[i][j] );
	}
    }

    return 10 * log10( num / den );
} /* image_snr */

/*--------------------------------------------------------------------------*/

/** \brief Calculates the Peak-Signal-to-Noise Ratio (PSNR) for a degraded  
           version with respect to its original image.

    Given an original image \f$g\f$ and a degraded version \f$f\f$, the 
    Peak-Signal-to-Noise Ratio (PSNR) is defined as
    \f[
      \mbox{PSNR}(f, g) := 10 \log_{10} 
      \left( 
        \frac{255^2}{\mbox{MSE}}
      \right) \enspace , \quad \mbox{and }
      \mbox{MSE} := 
      \frac{\sum_{i=1}^{nx} \sum_{j=1}^{ny} (f_{i,j} - g_{i,j})^2}
      {nx \cdot ny} \enspace .
    \f]

    Source: 
    C. Kervrann and J. Boulanger, Local adaptivity to variable smoothness
    for exemplar-based image regularization and representation.
    International Journal on Computer Vision, 2007, in print.
 */
double image_psnr(double **g,   /**< original image, unchanged  */
		  double **f,   /**< degraded image, unchagend */
		  long   nx,    /**< image size in x direction    */
		  long   ny,    /**< image size in y direction    */
		  long   bx,    /**< boundary size in x direction */
		  long   by)    /**< boundary size in y direction */
{
    long i, j;
    double den  = 0.0;

    for(i=bx; i<nx+bx; i++)
    {
	for(j=by; j<ny+by; j++)
	{
	    den = den + ( f[i][j] - g[i][j] ) * ( f[i][j] - g[i][j] );
	}
    }

    /* if original and degraded image coincide, we divide by zero! */
    return 10 * log10( 65025 * (double)(nx*ny) / den );
} /* image_snr */

/*--------------------------------------------------------------------------*/

/** \brief Calculates minimum, maximum, mean and variance of an image u.    
 */
void image_analyse(double **u,         /**< the image, unchanged         */
		   long   nx,          /**< image size in x direction    */
		   long   ny,          /**< image size in y direction    */
		   long   bx,          /**< boundary size in x direction */
		   long   by,          /**< boundary size in y direction */
		   double *min,        /**< minimum, output              */
		   double *max,        /**< maximum, output              */
		   double *mean,       /**< mean, output                 */ 
		   double *stdv)       /**< standard deviation, output   */
{
    long   i, j;       /* loop variable */
    double vari;       /* variance */
    double help;       /* auxiliary variable */
    double help2;      /* auxiliary variable */
    
    *min  = u[bx][by];
    *max  = u[bx][by];
    help2 = 0.0;

    for(i=bx; i<nx+bx; i++)
    {
	for(j=by; j<ny+by; j++)
	{
	    if(u[i][j] < *min) 
	    {
		*min = u[i][j];
	    }

	    if(u[i][j] > *max) 
	    {
		*max = u[i][j];
	    }

	    help2 += u[i][j];
	}
    }
    *mean = help2 / (double)(nx*ny);
    
    vari = 0.0;
    for(i=bx; i<nx+bx; i++)
    {
	for(j=by; j<ny+by; j++)
	{
	    help  = u[i][j] - *mean;
	    vari += help * help;
	}
    }
    vari  /= (double)(nx*ny);
    *stdv  = sqrt(vari);
    
    return;
} /* image_analyse */

/*--------------------------------------------------------------------------*/

/** \brief Calculates minimum, maximum, mean and variance of an image u and
           prints the result out on screen. 
 */
void image_analyse_print(double **u,    /**< input image, unchanged       */
			 long   nx,     /**< image size in x direction    */
			 long   ny,     /**< image size in y direction    */
			 long   bx,     /**< boundary size in x direction */
			 long   by,     /**< boundary size in y direction */
			 const char* name) /**< name of the image, only
					        used for screen output    */
{
    double min, max, mean, stdv;

    image_analyse(u, nx, ny, bx, by, &min, &max, &mean, &stdv);

    printf("analyse image %s: \n", name);
    printf("  min      : %11.6f\n", min);
    printf("  max      : %11.6f\n", max);
    printf("  mean     : %11.6f\n", mean);
    printf("  stdv     : %11.6f\n", stdv);
    printf("  l1norm   :   %13.6e\n", image_l1norm(u, nx, ny, bx, by) );
    printf("  l2norm   :   %13.6e\n", image_l2norm(u, nx, ny, bx, by) );
    printf("  linfnorm : %11.6f\n", image_linfnorm(u, nx, ny, bx, by) );

    return;
} /* image_analyse_print */

/*--------------------------------------------------------------------------*/

/** \brief Calculates minimum, maximum, mean and variance of an image u and
           prints the result out on screen. 
 */
void image_diff_analyse_print(
    double **g,    /**< original image, unchanged    */
    double **f,    /**< degraded image, unchanged    */
    long   nx,     /**< image size in x direction    */
    long   ny,     /**< image size in y direction    */
    long   bx,     /**< boundary size in x direction */
    long   by,     /**< boundary size in y direction */
    const char* name_g, /**< name of the original image, only
			     used for screen output           */
    const char* name_f) /**< name of the degraded image, only
			     used for screen output           */
{
    printf("analyse difference between \n");
    printf("  * original image %s \n", name_g);
    printf("  * degraded image %s \n", name_f);
    printf("  l1norm   :   %13.6e\n", 
	   image_diff_l1norm(g, f, nx, ny, bx, by) );
    printf("  l2norm   :   %13.6e\n", 
	   image_diff_l2norm(g, f, nx, ny, bx, by) );
    printf("  linfnorm : %11.6f\n", 
	   image_diff_linfnorm(g, f, nx, ny, bx, by) );
    printf("  SNR      : %11.6f\n", 
	   image_snr(g, f, nx, ny, bx, by) );

    return;
} /* image_diff_analyse_print */

/*--------------------------------------------------------------------------*/

#endif /* __IMANA_H__ */

