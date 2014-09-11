#include "watergrow.h"

/* -------------------------------------------------------------------------
Watershed transform for Image
Copyright (c) 2006 Oleg Sklyar
See: ../LICENSE for license, LGPL
------------------------------------------------------------------------- */

#include "tools.h"
#include <R_ext/Error.h>
#include <queue>
#include <map>

using namespace std;

class Pixel { 
public:
    double intensity;
    int x, y;
    double label;
    Pixel (double value, int _x, int _y, double l) : 
        intensity(value), x(_x), y(_y), label(l) {}
    bool underwater(double w,
                    double *kernel, PointXY * ksize,
                    double * labels, PointXY * dsize,
                    int min_votes ) {
        if (intensity > w) {
            return false;
        }
        else {
            if (min_votes == 0)
                return true;
            int votes = 0;
            int i, j, xx, yy, kcx, kcy;
            kcx = ksize->x / 2;
            kcy = ksize->y / 2;
            for ( i = -kcx; i <= kcx; i++ )
                for ( j = -kcy; j <= kcy; j++ ) {
                    if ( kernel[ (i + kcx) + (j + kcy) * ksize->x ] == 0) continue;
                    xx = x + i;
                    yy = y + j;
                    if ( xx < 0 || yy < 0 || xx >= dsize->x || yy >= dsize->y ) continue;
                    if ( labels[xx + yy * dsize->x] == label ) votes++;
                }
            return votes >= min_votes;
        }
    }
    void assign_nearest(double *kernel, PointXY * ksize,
                    double * labels, PointXY * dsize) {

        labels[ indexFromXY(x, y, dsize->x) ] = label;
        
        /*
        
        map<int,int> label_count;
        int neighbour;
        // initialize
        if (x > 0) {
            neighbour = indexFromXY(x - 1, y, dsize->x);
            if (labels[neighbour] > 0)
                label_count[ labels[neighbour] ] = 0;
        }
        if (y > 0) {
            neighbour = indexFromXY(x, y - 1, dsize->x);
            if (labels[neighbour] > 0)
                label_count[ labels[neighbour] ] = 0;
        }
        if (x < dsize->x - 1) {
            neighbour = indexFromXY(x + 1, y, dsize->x);
            if (labels[neighbour] > 0)
                label_count[ labels[neighbour] ] = 0;
        }
        if (y < dsize->y - 1) {
            neighbour = indexFromXY(x, y + 1, dsize->x);
            if (labels[neighbour] > 0)
                label_count[ labels[neighbour] ] = 0;
        }
        // count
        int i, j, xx, yy, kcx, kcy;
        kcx = ksize->x / 2;
        kcy = ksize->y / 2;
        for ( i = -kcx; i <= kcx; i++ )
            for ( j = -kcy; j <= kcy; j++ ) {
                if ( kernel[ (i + kcx) + (j + kcy) * ksize->x ] == 0) continue;
                xx = x + i;
                yy = y + j;
                if ( xx < 0 || yy < 0 || xx >= dsize->x || yy >= dsize->y ) continue;
                if ( label_count.count( labels[xx + yy * dsize->x] ) > 0 )
                    label_count[ labels[xx + yy * dsize->x] ]++;
            }
        // get maximum
        map<int,int>::iterator it;
        int label_maximum = 0, maximum = 0;
        for ( it=label_count.begin() ; it != label_count.end(); it++ ) {
            if ((*it).second > maximum ) {
                label_maximum = (*it).first;
                maximum = (*it).second;
            }
        }
        // assign
        labels[ indexFromXY(x, y, dsize->x) ] = label_maximum;
*/
    }
};

bool
add_candidates( queue<Pixel> &pq,
                double min_water, double max_water,
                double *image,
                int x, int y, int nx, int ny,
                double label,
                double *labels ) {

    int neighbour;
    bool changed = false;
    
    /* 4-connected */
    if (x > 0) {
        neighbour = indexFromXY(x - 1, y, nx);
        if ( (labels[neighbour] == 0) // not labelled nor seen
             && (image[neighbour] <= max_water)
             && (image[neighbour] >= min_water) ) {
            pq.push(Pixel(image[neighbour], x - 1, y, label));
            labels[neighbour] = -1; // seen
            changed = true;
        }
    }
    if (y > 0) {
        neighbour = indexFromXY(x, y - 1, nx);
        if ( (labels[neighbour] == 0) // not labelled nor seen
             && (image[neighbour] <= max_water)
             && (image[neighbour] >= min_water) ) {
            pq.push(Pixel(image[neighbour], x, y - 1, label));
            labels[neighbour] = -1; // seen
            changed = true;
        }
    }
    if (x < (nx-1)) {
        neighbour = indexFromXY(x + 1, y, nx);
        if ( (labels[neighbour] == 0) // not labelled nor seen
             && (image[neighbour] <= max_water)
             && (image[neighbour] >= min_water) ) {
            pq.push(Pixel(image[neighbour], x + 1, y, label));
            labels[neighbour] = -1; // seen
            changed = true;
        }
    }
    if (y < (ny-1)) {
        neighbour = indexFromXY(x, y + 1, nx);
        if ( (labels[neighbour] == 0) // not labelled nor seen
             && (image[neighbour] <= max_water)
             && (image[neighbour] >= min_water) ) {
            pq.push(Pixel(image[neighbour], x, y + 1, label));
            labels[neighbour] = -1; // seen
            changed = true;
        }
    }
    return changed;
}

/*----------------------------------------------------------------------- */
SEXP
watergrow (SEXP _map, SEXP _labels, SEXP _min_water, SEXP _max_water,
           SEXP _step_size, SEXP _kernel, SEXP _min_votes) {
    SEXP res;
    int i, nx, ny, count,min_votes, nprotect = 0;
    double *kernel;
    double min_water, max_water, step_size;
    PointXY pt, size, ksize;

    nx = size.x = INTEGER ( GET_DIM(_map) )[0];
    ny = size.y =  INTEGER ( GET_DIM(_map) )[1];
    min_water = REAL( _min_water )[0];
    max_water = REAL( _max_water )[0];
    step_size = REAL( _step_size )[0];
    min_votes = INTEGER ( _min_votes )[0];

    kernel = &( REAL(_kernel)[ 0 ] );
    ksize.x = INTEGER ( GET_DIM(_kernel) )[0];
    ksize.y = INTEGER ( GET_DIM(_kernel) )[1];

    PROTECT ( res = Rf_duplicate(_labels) );
    nprotect++;
    
    double * map = &( REAL(_map)[  0 ] );
    double * labels = &( REAL(res)[  0 ] );

    queue<Pixel> candidates;

    for ( i = 0; i < nx * ny; i++ ) {
        if ( labels[ i ] > 0 ) { // 0 for background, -1 for seen
            pt = pointFromIndex( i, nx );
            add_candidates( candidates,
                            min_water, max_water,
                            map,
                            pt.x, pt.y,
                            nx, ny,
                            labels[i], labels );
        }
    }

    double water_level = min_water;
    bool changed;
    while ( water_level <= max_water + step_size ) {

        R_CheckUserInterrupt();
        changed = true;

        while( changed ) {
            changed = false;
            count = candidates.size();
            while (count > 0 ) {
                count--;
                Pixel current = candidates.front();
                candidates.pop();
                if (current.underwater(water_level,
                                       kernel, &ksize,
                                       labels, &size,
                                       min_votes)) {
                    i = indexFromXY(current.x, current.y, nx);
                    labels[i] = current.label;
                    changed = add_candidates( candidates,
                                              min_water, max_water,
                                              map,
                                              current.x, current.y,
                                              nx, ny,
                                              current.label, labels )
                        || changed;
                }
                else {
                    candidates.push(current);
                }
            }
        } // while something changes at the same water level
        
        if (water_level <= max_water + step_size) {
            water_level += step_size;
        }
    }

    while (! candidates.empty() ) {
        Pixel current = candidates.front();
        candidates.pop();
        i = indexFromXY(current.x, current.y, nx);
        labels[i] = current.label;            
        add_candidates( candidates,
                        min_water, max_water,
                        map,
                        current.x, current.y,
                        nx, ny,
                        current.label, labels );
    }

    UNPROTECT (nprotect);
    return res;
}
