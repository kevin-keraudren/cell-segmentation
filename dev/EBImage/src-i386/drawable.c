#include "drawable.h"

/* -------------------------------------------------------------------------
Drawing on images
Copyright (c) 2007 Oleg Sklyar
See: ../LICENSE for license, LGPL
------------------------------------------------------------------------- */

#include "tools.h"
#include "conversions.h"

#include <R_ext/Error.h>
/* this module uses MagickWand API instead of MagickCore */
#include <wand/magick-wand.h>

SEXP
lib_drawText (SEXP obj, SEXP xylist, SEXP textlist, SEXP thefont, SEXP thecol) {
  SEXP res, xy, text;
  int i, im, nz, mode, nval, nappended=0, nprotect=0, colcount=0;
  DrawingWand * dwand;
  MagickWand  * mwand;
  PixelWand * pwand;
  double * dxy;
  Image * images, * image, * newimages;
  const char * str;
  
  validImage(obj,0);

  mode = getColorMode(obj);
  images = sexp2Magick(obj);
  nz = GetImageListLength(images);
  
  if ( LENGTH(xylist) != LENGTH(textlist) || LENGTH(xylist) != nz )
    error("lists of coordinates 'xy' labels 'labels' must be of the same length as the number of frames");

  newimages = NewImageList();

  /* create empty wand */
  dwand = NewDrawingWand();
  pwand = NewPixelWand();
  /* loop through images */
  for ( im = 0; im < nz; im++ ) {
    /* create magick wand from one image, this does NOT copy */
    mwand = NewMagickWandFromImage(GetFirstImageInList(images));
    /* get pointers to current list elements if those are lists */
    xy = VECTOR_ELT(xylist, im);
    dxy = REAL(xy);
    text = VECTOR_ELT(textlist, im);
    nval = LENGTH(text);
    if ( nval > 0 && LENGTH(xy) >= 2 * nval ) {
      /* clear drawing wand */
      ClearDrawingWand(dwand);
      str = CHAR(STRING_ELT(VECTOR_ELT(thefont, 0),0));
      if (str)
        DrawSetFontFamily(dwand, str);
      switch ( INTEGER(VECTOR_ELT(thefont, 1))[0]) {
        case 1:  DrawSetFontStyle(dwand, ItalicStyle); break;
        case 2:  DrawSetFontStyle(dwand, ObliqueStyle); break;
        default: DrawSetFontStyle(dwand, NormalStyle);
      }
      DrawSetStrokeAntialias(dwand, INTEGER(VECTOR_ELT(thefont, 4))[0]);
      DrawSetFontWeight(dwand, REAL(VECTOR_ELT(thefont, 3))[0]);
      DrawSetFontSize(dwand, REAL(VECTOR_ELT(thefont, 2))[0]);
      PixelSetColor(pwand, CHAR(STRING_ELT(thecol, colcount)));
      colcount++;
      if ( colcount >= LENGTH(thecol) ) colcount = 0;
      DrawSetFillColor(dwand, pwand);
      /* no need: DrawSetStrokeColor(dwand, PixelWand); */
      /* add text to the wand */
      for ( i = 0; i < nval; i++ ) {
        str = CHAR(STRING_ELT(text, i));
        if (str)
          DrawAnnotation(dwand, dxy[i], dxy[i + nval], (unsigned char *)str );
      }
      /* draw the wand */
      MagickDrawImage(mwand, dwand);
    }
    else {
      /* do not draw if more text labels than coordinates */
      warning("not enough coordinate points to output all labels");
    }
    image = GetImageFromMagickWand(mwand);
    AppendImageToList (&newimages, image);
    if ( nappended == 0 ) {
      /* copy attributes once and only if image contains more than 1 element */
      nappended = 1;
      newimages->compression = images->compression;
      strcpy (newimages->filename, image->filename);
      newimages->x_resolution = images->x_resolution;
      newimages->y_resolution = images->y_resolution;
    }
    /* the image is destroyed below in DestroyImage 
    mwand = DestroyMagickWand(mwand);
    */
    image = GetFirstImageInList(images);
    RemoveFirstImageFromList(&images);
    image = DestroyImage(image);
  }
  pwand = DestroyPixelWand(pwand);
  dwand = DestroyDrawingWand(dwand);
  
  images = DestroyImageList(images);

  PROTECT( res = magick2SEXP(newimages, mode) );
  nprotect++;
  newimages = DestroyImageList(newimages);

  UNPROTECT( nprotect );
  return res;
}

#define SET_PIXEL(a, width, height, x, y, color) if (((x)>=0) & ((x)<(width)) & ((y)>=0) & ((y)<(height))) a[(x) + (y)*(width)] = (color);
void rasterCircle(double *a, int width, int height, int x0, int y0, int radius, double color, int fill) {
  int f = 1 - radius;
  int ddF_x = 1;
  int ddF_y = -2 * radius;
  int x = 0;
  int y = radius;
  int i;
  
  if (fill) {
    for (i=x0-radius; i<=x0+radius ;i++) SET_PIXEL(a, width, height, i, y0, color);
    for (i=y0-radius; i<=y0+radius ;i++) SET_PIXEL(a, width, height, x0, i, color);
  } else {
    SET_PIXEL(a, width, height, x0, y0 + radius, color);
    SET_PIXEL(a, width, height, x0, y0 - radius, color);
    SET_PIXEL(a, width, height, x0 + radius, y0, color);
    SET_PIXEL(a, width, height, x0 - radius, y0, color);
  }
  
  while(x < y) {
    if(f >= 0) {
      y--;
      ddF_y += 2;
      f += ddF_y;
    }
    x++;
    ddF_x += 2;
    f += ddF_x;
    if (fill) {
      for (i=x0-x; i<=x0+x ;i++) SET_PIXEL(a, width, height, i, y0+y, color);
      for (i=x0-x; i<=x0+x ;i++) SET_PIXEL(a, width, height, i, y0-y, color);
      for (i=x0-y; i<=x0+y ;i++) SET_PIXEL(a, width, height, i, y0+x, color);
      for (i=x0-y; i<=x0+y ;i++) SET_PIXEL(a, width, height, i, y0-x, color);
    }
    else {
      SET_PIXEL(a, width, height, x0 + x, y0 + y, color);
      SET_PIXEL(a, width, height, x0 - x, y0 + y, color);
      SET_PIXEL(a, width, height, x0 + x, y0 - y, color);
      SET_PIXEL(a, width, height, x0 - x, y0 - y, color);
      SET_PIXEL(a, width, height, x0 + y, y0 + x, color);
      SET_PIXEL(a, width, height, x0 - y, y0 + x, color);
      SET_PIXEL(a, width, height, x0 + y, y0 - x, color);
      SET_PIXEL(a, width, height, x0 - y, y0 - x, color);
    }
  }
}

// draw a circle on the 2D image _a using (x, y, z, radius) and color (red, green, blue)
// if colormode = Grayscale, only the red component is used
SEXP drawCircle(SEXP _a, SEXP _xyzr, SEXP _rgb, SEXP _fill) {
  SEXP _res;
  int nprotect = 0;
  int width, height;
  int x, y, z, radius;
  int redstride, greenstride, bluestride;
  double *res;
  int fill;

  // check image validity and copy _a
  validImage(_a, 0);
  PROTECT(_res=Rf_duplicate(_a));
  nprotect++;

  width = INTEGER(GET_DIM(_res))[0];
  height = INTEGER(GET_DIM(_res))[1];

  // get strides
  x = INTEGER(_xyzr)[0];
  y = INTEGER(_xyzr)[1];
  z = INTEGER(_xyzr)[2];
  radius = INTEGER(_xyzr)[3];
  fill = INTEGER(_fill)[0];
  getColorStrides(_res, z, &redstride, &greenstride, &bluestride);
  res = REAL(_res);

  // draw circle
  if (getColorMode(_res)==MODE_GRAYSCALE) {
    rasterCircle(&res[redstride], width, height, x, y, radius, REAL(_rgb)[0], fill);
  } 
  else if (getColorMode(_res)==MODE_COLOR) {
    rasterCircle(&res[redstride], width, height, x, y, radius, REAL(_rgb)[0], fill);
    rasterCircle(&res[greenstride], width, height, x, y, radius, REAL(_rgb)[1], fill);
    rasterCircle(&res[bluestride], width, height, x, y, radius, REAL(_rgb)[2], fill);
  }
  
  UNPROTECT (nprotect);
  return _res;
}

// taken from biOps package
SEXP
imageConvolve (SEXP obj, SEXP _kern) {
    SEXP res;
    int x, y, z, i, j, nx,ny,nz,kx,ky,kcx,kcy,xx,yy;
    double * src, *tgt, *kernel, value;

    validImage(obj,0);

    nx = INTEGER ( GET_DIM(obj) )[0];
    ny = INTEGER ( GET_DIM(obj) )[1];
    nz = getNumberOfFrames(obj,0);

    kernel = REAL (_kern);
    kx = INTEGER ( GET_DIM(_kern) )[0];
    ky = INTEGER ( GET_DIM(_kern) )[1];

    kcx = kx / 2;
    kcy = ky / 2;

    int nprotect = 0;

    PROTECT ( res = Rf_duplicate(obj) );
    nprotect++;

    for ( z = 0; z < nz; z++ ) {
        tgt = &( REAL(res)[z * nx * ny] );
        src = &( REAL(obj)[z * nx * ny] );

        for (x = 0; x < nx; x++){
            for (y = 0; y < ny; y++){
                value = 0;
                for (i = -kcx; i <= kcx; i++){
                    for (j = -kcy; j <= kcy; j++){
                        xx = x + i;
                        yy = y + j;
                        if ( xx < 0 || yy < 0 || xx >= nx || yy >= ny )continue;
                        value += src[ xx + yy*nx] * kernel[ (i + kcx) + (j + kcy) * kx];
                    }
                }
                tgt[x + y*nx] = value;
            }
        }
    }
    UNPROTECT (nprotect);    
    return res;
}

SEXP
imageReplace (SEXP obj, SEXP _colors, SEXP _do_zero) {
    SEXP res;
    int z, j, nx,ny,nz,do_zero;
    double *tgt, *colors;
    double *src;

    validImage(obj,0);

    nx = INTEGER ( GET_DIM(obj) )[0];
    ny = INTEGER ( GET_DIM(obj) )[1];
    nz = getNumberOfFrames(obj,0);

    colors = REAL (_colors);
    do_zero = INTEGER ( _do_zero)[0];
    //nb_colors = INTEGER ( GET_DIM(_colors) )[0];

    int nprotect = 0;

    PROTECT ( res = Rf_duplicate(obj) );
    nprotect++;

    for ( z = 0; z < nz; z++ ) {
        tgt = &( REAL(res)[z * nx * ny] );
        src = &( REAL(obj)[z * nx * ny] );

        if (do_zero) {
            for (j=0; j< nx*ny; j++)
                tgt[ j ] = colors[ (int)(src[j]) ];
        }
        else {
            for (j=0; j< nx*ny; j++) {
                if ( src[j] > 0 )
                    tgt[ j ] = colors[ (int)(src[j]) - 1 ];
            }
        }

    }
    
    UNPROTECT (nprotect);    
    return res;
}

SEXP
lib_drawPolyline (SEXP obj, SEXP xylist, SEXP stroke_color, SEXP stroke_width,
                  SEXP fill_color, SEXP antialias ) {
  SEXP res, xy, text;
  int i, im, nz,nxy, mode, nval, nappended=0, nprotect=0, colcount=0,npoly,p;
  DrawingWand * dwand;
  MagickWand  * mwand;
  PixelWand * pwand_stroke, * pwand_fill;
  double * dxy;
  Image * images, * image, * newimages;
  const char * str;
  
  validImage(obj,0);

  mode = getColorMode(obj);
  images = sexp2Magick(obj);
  nz = GetImageListLength(images);
  npoly = LENGTH(xylist);

  newimages = NewImageList();

  /* create empty wand */
  dwand = NewDrawingWand();
  pwand_stroke = NewPixelWand();
  pwand_fill = NewPixelWand();
  /* loop through images */
  for ( im = 0; im < nz; im++ ) {
    /* create magick wand from one image, this does NOT copy */
    mwand = NewMagickWandFromImage(GetFirstImageInList(images));
    /* get pointers to current list elements if those are lists */
    for (p=0;p<npoly;p++) {
    xy = VECTOR_ELT(xylist, p);
    nxy = LENGTH(xy) / 2;
    dxy = REAL(xy);
    /* nxy=INTEGER(GET_DIM(xylist))[0]; */
    /* dxy = REAL(xylist); */
      /* clear drawing wand */
      ClearDrawingWand(dwand);
	  ClearPixelWand(pwand_stroke);
	  ClearPixelWand(pwand_fill);

      DrawSetStrokeAntialias(dwand, MagickFalse);
      DrawSetStrokeWidth(dwand,REAL(stroke_width)[0]);
	  DrawSetFillOpacity(dwand, 0.0);
	  DrawSetStrokeOpacity(dwand, 1.0);

      PixelSetColor(pwand_stroke, CHAR(STRING_ELT(stroke_color, colcount)));
	  PixelSetOpacity(pwand_stroke, 1.0);
      PixelSetColor(pwand_fill, "none");
	  PixelSetOpacity(pwand_fill, 0.0);
      //PixelSetColor(pwand_fill, CHAR(STRING_ELT(fill_color, colcount)));
       
      /* colcount++; */
      /* if ( colcount >= LENGTH(thecol) ) colcount = 0; */
      
      DrawSetStrokeColor(dwand, pwand_stroke);
      //DrawSetFillColor(dwand, pwand_fill);
      DrawPolyline(dwand, nxy, dxy );
      /* draw the wand */
      MagickDrawImage(mwand, dwand);
    }

    image = GetImageFromMagickWand(mwand);
    AppendImageToList (&newimages, image);
    if ( nappended == 0 ) {
      /* copy attributes once and only if image contains more than 1 element */
      nappended = 1;
      newimages->compression = images->compression;
      strcpy (newimages->filename, image->filename);
      newimages->x_resolution = images->x_resolution;
      newimages->y_resolution = images->y_resolution;
    }
    /* the image is destroyed below in DestroyImage 
    mwand = DestroyMagickWand(mwand);
    */
    image = GetFirstImageInList(images);
    RemoveFirstImageFromList(&images);
    image = DestroyImage(image);    
  }
  pwand_stroke = DestroyPixelWand(pwand_stroke);
  pwand_fill = DestroyPixelWand(pwand_fill);
  dwand = DestroyDrawingWand(dwand);
  
  images = DestroyImageList(images);

  PROTECT( res = magick2SEXP(newimages, mode) );
  nprotect++;
  newimages = DestroyImageList(newimages);

  UNPROTECT( nprotect );
  return res;
}
