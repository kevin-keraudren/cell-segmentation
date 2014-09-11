#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include "ocontour.h"

#define MAX_NB_POINTS 32768

static int rotr[8]={-1,-1,-1, 0, 1, 1, 1, 0};
static int rotc[8]={-1, 0, 1, 1, 1, 0,-1,-1};
static int dir [9]={5,4,3,6,-1,2,7,0,1};

/* dir */
/*  0 */
/* 3 1 */
/*  2 */
int on_right(int i,int j,int dir, int *data, int k, int nx) {
    switch (dir) {
    case 0: return data[i+1 + nx*j] == k; break;
    case 1: return data[i+nx*(j+1)] == k; break;
    case 2: return data[i-1+nx*j] == k; break;
    case 3: return data[i+nx*(j-1)] == k; break;
    }
}
int in_front(int i,int j,int dir, int *data, int k, int nx) {
    switch (dir) {
    case 0: return data[i + nx*(j-1)] == k; break;
    case 1: return data[i+1+nx*j] == k; break;
    case 2: return data[i+nx*(j+1)] == k; break;
    case 3: return data[i-1+nx*j] == k; break;
    }
}

SEXP ocontour(SEXP _image, SEXP _external) {
    int *image, *padded_image = NULL;
    int i, j, k, direction, nbCells, external;
    int x,y , x0, y0, ndirection, x2,y2, nx,ny;
    int nprotect=0, nboc,  *octemp,moved;
    SEXP _res, _oc;
  
    // Transfer variables
    nx=INTEGER(GET_DIM(_image))[0];
    ny=INTEGER(GET_DIM(_image))[1];
    image=INTEGER(_image);
    external=INTEGER(_external)[0];
  
    // Compute number of objects
    nbCells=0;
    for (i=0; i<nx*ny; i++) {
        if (image[i]>nbCells) nbCells=image[i];
    }
    nbCells++;

    if (external) {
        nx +=2;
        ny += 2;
        padded_image =(int *)R_Calloc(nx*ny, int);
        for ( i=0; i < nx; i++) {
            for (j=0; j< ny; j++) {
                if (i==0||j==0||i == nx -1||j==ny-1)
                    padded_image[i+nx*j] = 0;
                else
                    padded_image[i+nx*j] = image[i-1+(nx-2)*(j-1)];
            }
        }
        image = padded_image;
    }
 
    // Output result
    PROTECT(_res = allocVector(VECSXP, nbCells));
    nprotect++;

    // Temporary vector to store the current oriented contour
    octemp=(int *)R_Calloc(MAX_NB_POINTS*2, int);

    // For each object, except the 0-th one (background)
    for (k=1; k<nbCells; k++) {
        nboc=0;

        // Find min (r,c) starting point for object k
        i=0;
        while (image[i]!=k && i<nx*ny) i++;
        if (i!=nx*ny) {
            x=i%nx;
            y=i/nx;

            // Starting points of the oriented contour
            if (external) {
                x0=x;
                y0=y=y-1;
                direction=1;
            }
            else {
                x0=x;
                y0=y;
            }

            // Turn around the edges of the object
            if (! external)
                direction=0;
            do {
                if (nboc % 1000)
                    R_CheckUserInterrupt();
                // Stores (r,c) in the oriented contour matrix
                // R indices start from 1
                if (external) {
                    if ( x == 0 )
                        octemp[2*nboc] = 1;
                    else if (x == nx - 1)
                        octemp[2*nboc] = nx -2;
                    else
                        octemp[2*nboc] = x;
                    if (y == 0)
                        octemp[2*nboc+1] = 1;
                    else if (y== ny -1)
                        octemp[2*nboc+1] = ny - 2;
                    else
                        octemp[2*nboc+1] = y;           
                }
                else {
                    octemp[2*nboc] = x + 1;
                    octemp[2*nboc+1] = y + 1;
                }
                if (nboc<MAX_NB_POINTS) nboc++;
	
                // Change direction
                if (!external) {
                    for (j=0; j<8; j++) {
                        ndirection=(j+direction)%8;
                        x2=x+rotr[ndirection];
                        y2=y+rotc[ndirection];
                        if (x2>=0 && y2>=0 && x2<nx && y2<ny) {
                            if (image[x2+y2*nx]==k)  break;
                        }
                    }
                    if (j!=8) {
                        direction=dir[(x2-x+1)+3*(y2-y+1)];
                        x=x2;
                        y=y2;
                    }            
                }
                else {
                    moved = 0;
                    while ((! moved)) {
                        if (on_right(x,y,direction, image, k, nx)) {
                            if (in_front(x,y,direction, image, k, nx)) {
                                // turn left
                                // sign issue when using the % (modulo)
                                if (direction == 0)
                                    direction = 3;
                                else
                                    direction--;
                            }
                            else {
                                // go straight
                                switch (direction) {
                                case 0: y -= 1; break;
                                case 1: x += 1; break;
                                case 2: y += 1; break;
                                case 3: x -= 1; break;
                                }
                                moved = 1;
                            }
                        }
                        else {
                            // turn right
                            switch (direction) {
                            case 0: x += 1; break;
                            case 1: y += 1; break;
                            case 2: x -= 1; break;
                            case 3: y -= 1; break;
                            }
                            direction = (direction +1) % 4;
                            moved = 1;
                        }
                    }
                }
            } while (x!=x0 || y!=y0);
        }
    
        // Copy octemp in an element of _res
        _oc=allocVector(INTSXP, nboc*2);
        SET_VECTOR_ELT(_res, k, _oc);
        memcpy(INTEGER(_oc), octemp, nboc*2*sizeof(int));
    } // k

    // Free oct
    R_Free(octemp);
    if (external) {
       image = NULL;
       R_Free(padded_image);
    }

    UNPROTECT (nprotect);
    return(_res);
}
