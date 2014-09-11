# Drawables for class Image

# Copyright (c) 2007 Oleg Sklyar

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.

# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

# See the GNU Lesser General Public License for more details.
# LGPL license wording: http://www.gnu.org/licenses/lgpl.html

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
drawfont = function(family=switch(.Platform$OS.type, windows="Arial", "helvetica"),
  style="n", size=14, weight=200, antialias=TRUE) {
  res <- list(family=family, style=style, size=size, weight=weight, antialias=antialias)
  class(res) <- "DrawFont"
  res
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
drawtext = function(img, xy, labels, font, col) {   
  validImage(img)
  if (is.numeric(xy)) {
    xy=list(as.numeric(xy))
    labels=list(labels)
  }
  if (missing(font)) font=drawfont()
  if (missing(col)) col="white"
  
  ## Mac stop: imageMagick/MagickWand/drawtext doesn't understand MacOS fonts
  if (length(grep("apple", Sys.getenv("R_PLATFORM")))>0) {
    warning('drawtext doesn\'t work on Mac OS due to ImageMagick/Mac OS fonts incompatibility')
    return(img)
  }
    
  if (length(xy) != length(labels) || length(xy) != getNumberOfFrames(img,'render'))
    stop("lists of coordinates 'xy' labels 'labels' must be of the same length as the number of render frames")
  xy <- lapply(xy, as.numeric)
  for ( i in seq_along(labels)) 
    if (!is.character(labels[[i]]))
      stop("all elements of 'labels' must be of class 'character'")
  if ( !is(font, "DrawFont") )
      stop("to set the font use the 'drawfont' function which returns an S3 class 'DrawFont', modify the slots as needed")

  font$style <- as.integer(switch(tolower(substr(font$style,1,1)), i=1, o=2, 0))
  font$size <- as.numeric(font$size)
  font$weight <- as.numeric(font$weight)
  font$antialias <- as.logical(font$antialias)
  return(.Call("lib_drawText", castImage(img), xy, labels, font, col, PACKAGE='EBImage'))
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
drawPolyline = function(img, xy, stroke.color="white",stroke.width=1,fill.color="none",antialias=FALSE) {   
  validImage(img)
  
  if (is.numeric(xy)) {
    xy=list(as.numeric(xy))
  }

  # FIXME: the farfetched code below does:
  # - substract 1 to coordinates to translate from R coordinates to C
  # coordinates
  # - returns everything 'as.numeric'
  # - keeping the structure of a list of lists...
  for (i in 1:length(xy)) {
    xy[[i]] <- do.call('rbind', lapply(1:dim(xy[[i]])[1], function(j){
     rbind(as.numeric(xy[[i]][j,1] - 1),as.numeric(xy[[i]][j,2] - 1) )
 }))
  }

  return(.Call("lib_drawPolyline", castImage(img), xy,
               stroke.color, stroke.width, fill.color, as.integer(antialias), PACKAGE='EBImage'))
}

drawCircle = function(img, x, y, radius, col, fill=FALSE, z=1) {
  validImage(img)
  if (colorMode(img)==TrueColor) stop("this method doesn't support the \'TrueColor\' color mode")
  if (any(is.na(img))) stop("'x' shouldn't contain any NAs")

  ## check whether parameters are OK
  if (missing(radius)) stop("'radius' is missing")
  if (radius<1) stop("'radius' must be positive integer")
  if (z<1 | z>getNumberOfFrames(img, 'render')) stop("'z' must be a positive integer lower than the number of image frames")
  xyzr = as.integer(c(x, y, z-1, radius))
  if (length(xyzr)!=4 || any(is.na(xyzr))) stop("'x', 'y', 'z' and 'radius' must be scalar values")
  fill = as.integer(fill)
  if (length(fill)!=1)  stop("'fill' must be a logical")
  
  if (colorMode(img)==Color) {
    rgb = as.numeric(col2rgb(col)/255)
    if (length(rgb)!=3 || any(is.na(rgb))) stop("In Color mode, 'col' must be a valid color")
  } else {
    rgb = as.numeric(c(col, 0, 0))
    if (length(rgb)!=3 || any(is.na(rgb))) stop("In Grayscale mode, 'col' must be a scalar value")
  }
  
  invisible(.Call("drawCircle", castImage(img), xyzr, rgb, fill, PACKAGE='EBImage'))
}

imageReplace = function(x,colors, do.zero=FALSE) {

  if ( max(x) > length(colors)) stop("not enough colors")

  .Call("imageReplace", castImage(x), as.numeric(colors), as.integer(do.zero), PACKAGE='EBImage')
}

colorize = function(img,colormap="jet") {
  if ( ! is.character(colormap)) stop("colormap must be a 'character'")

  d = dim(img)
  
  if (colormap == "jet") {
    jet.colors <-
      colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                         "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    palette = jet.colors(256)
  }
  else if (colormap == "rainbow") {
    palette = rainbow(256)
  }
  else stop("Unknown colormap")

  datac = col2rgb(palette)/255
  img = normalize(img)*255
  rgbImage(imageReplace(img, datac[1,]),
imageReplace(img, datac[2,]),imageReplace(img, datac[3,]) )
  
}
