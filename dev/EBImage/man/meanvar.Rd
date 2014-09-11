\name{meanvar}

\alias{meanvar}

\title{Mean intensity and variance of objects}

\description{
  Computes the mean intensity and variance of image objects.
}

\usage{
  meanvar(x, ref, nobj = max(x)) 
}

\arguments{
  \item{x}{An \code{Image} object or an array containing object masks.
    Object masks are sets of pixels with the same unique integer value.}
    
  \item{ref}{An \code{Image} object or an array, containing the
    intensity values of the objects.}
  
  \item{nobj}{The total number of objects: in the case of subobjects,
  there can be none corresponding to a superobject. Giving a value for
  \code{nobj} assures you that the resulting list will have the correct length.}

}

\value{
  A matrix of features 'mean', 'var' and 'size' (in pixels).
}

\seealso{
  \code{\link{getFeatures}}, \code{\link{bwlabel}}, \code{\link{watershed}}, \code{\link{propagate}}
}

\author{
  Kevin Keraudren, \email{kevin.keraudren10@imperial.ac.uk}
}

\examples{
  ## load cell nucleus images
  x = readImage(system.file('images', 'nuclei.tif', package='EBImage'))
  if (interactive()) display(x)

  ## computes object mask
  y = thresh(x, 10, 10, 0.05)
  y = opening(y, makeBrush(5, shape='disc'))
  mask = fillHull(bwlabel(y))
  if (interactive()) display(mask, title='Cell nuclei')

  ## mean and variance
  mvar = meanvar(mask[,,1], x[,,1])
}


