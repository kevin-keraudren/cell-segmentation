#include <R_ext/Error.h>
#include "tools.h"
#include "morphology.h"

/* -------------------------------------------------------------------------
Morphological filters for Image
Copyright (c) 2006 Oleg Sklyar
See: ../LICENSE for license, LGPL
------------------------------------------------------------------------- */

#define ERODE  0
#define DILATE 1

int _match (numeric * kernel, PointXY * ksize, numeric * data, PointXY * dsize,
            PointXY * at, numeric mismatch);

numeric
_greymatch (numeric * kernel, PointXY * ksize, numeric * data, PointXY * dsize,
            PointXY * at, int what, numeric _min, numeric _max);

/*----------------------------------------------------------------------- */
SEXP
lib_erode_dilate (SEXP x, SEXP kernel, SEXP what, SEXP binary) {
    numeric resetTo, * tgt, * src, *kern, min, max;
    int nz, i, j, nprotect;
    int * dim;
    PointXY size, ksize, pt;
    SEXP res;

    validImage(x,0);
    validImage(kernel,0);

    /* value to reset the checked part t */
    if ( INTEGER(what)[0] == DILATE )
        resetTo = 1.0; /* checking background, reseting to 1 */
    else
        resetTo = 0.0; /* checking foreground, reseting to 0 */
    dim = INTEGER ( GET_DIM(x) );
    size.x = dim[0];
    size.y = dim[1];
    nz = getNumberOfFrames(x,0);

    kern = REAL (kernel);
    ksize.x = INTEGER ( GET_DIM(kernel) )[0];
    ksize.y = INTEGER ( GET_DIM(kernel) )[1];
    nprotect = 0;

    PROTECT ( res = Rf_duplicate(x) );
    nprotect++;

    for ( i = 0; i < nz; i++ ) {
        tgt = &( REAL(res)[i * size.x * size.y] );
        src = &( REAL(x)[i * size.x * size.y] );

        if ( ! INTEGER(binary)[0] ) {
            min = max = src[0];
            for ( j = 0; j < size.x * size.y; j++ ) {
                if (src[j] > max)
                    max = src[j];
                if (src[j] < min)
                    min = src[j];
            }
            for ( j = 0; j < size.x * size.y; j++ ) {
                pt = pointFromIndex (j, size.x);
                tgt[j] = _greymatch(kern, &ksize, src, &size, &pt,
                                    INTEGER(what)[0], min , max);
            }
        }
        else {
            for ( j = 0; j < size.x * size.y; j++ ) {
                if ( tgt[j] == resetTo ) continue;
                pt = pointFromIndex (j, size.x);
                if ( !_match(kern, &ksize, src, &size, &pt, resetTo) )
                    tgt[j] = resetTo;
            }
        }
    }
    
    UNPROTECT (nprotect);
    return res;
}


/*----------------------------------------------------------------------- */
int
_match (numeric * kernel, PointXY * ksize, numeric * data, PointXY * dsize, PointXY * at, numeric mismatch) {
    int i, j, xx, yy, kcx, kcy;

    kcx = ksize->x / 2;
    kcy = ksize->y / 2;
    for ( i = -kcx; i <= kcx; i++ )
        for ( j = -kcy; j <= kcy; j++ ) {
            if ( kernel[ (i + kcx) + (j + kcy) * ksize->x ] == 0) continue;
            xx = at->x + i;
            yy = at->y + j;
            if ( xx < 0 || yy < 0 || xx >= dsize->x || yy >= dsize->y ) continue;
            if ( data[xx + yy * dsize->x] == mismatch ) return 0;
        }
    return 1;
}

/*----------------------------------------------------------------------- */
numeric
_greymatch (numeric * kernel, PointXY * ksize, numeric * data, PointXY * dsize,
            PointXY * at, int what, numeric _min, numeric _max) {
    int i, j, xx, yy, kcx, kcy;
    numeric value;

    numeric min = _max;
    numeric max = _min;

    kcx = ksize->x / 2;
    kcy = ksize->y / 2;
    for ( i = -kcx; i <= kcx; i++ )
        for ( j = -kcy; j <= kcy; j++ ) {
            if ( kernel[ (i + kcx) + (j + kcy) * ksize->x ] == 0) continue;
            xx = at->x + i;
            yy = at->y + j;
            if ( xx < 0 || yy < 0 || xx >= dsize->x || yy >= dsize->y )continue;
            if (what == ERODE) {
                value = data[xx + yy * dsize->x] - kernel[ (i + kcx) + (j + kcy)
    * ksize->x ];
                if (value < min)
                    min = value;
            }
            else {
                value = data[xx + yy * dsize->x] + kernel[ (i + kcx) + (j + kcy)
    * ksize->x ];
                if (value > max)
                    max = value;
            }
        }
    if (what == ERODE)
        return min;
    else
        return max;
}

