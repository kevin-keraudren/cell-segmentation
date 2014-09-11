\name{ocontour}

\alias{ocontour}
\title{Oriented contours}

\description{
  Computes the oriented contour of objects.
}

\usage{
  ocontour(x, external=FALSE)
}

\arguments{
  \item{x}{An \code{Image} object or an array, containing objects.
    Only integer values are considered. Pixels of value 0 constitute
    the background. Each object is a set of pixels with the same unique
    integer value. Objects are assumed connected.}

  \item{external}{A logical indicating whether to return the internal or
  external contour.}

}

\value{
  A list of matrices, containing the coordinates of object oriented contours.
}

\author{
  Gregoire Pau, \email{gpau@ebi.ac.uk}, 2008
}

\examples{
  x = readImage(system.file("images", "shapes.png", package="EBImage"))
  x = x[1:120,50:120]
  if(interactive()) display(x)
  oc = ocontour(x)
  plot(oc[[1]], type='l')
  points(oc[[1]], col=2)
  oc2 = ocontour(x,external=TRUE)
  lines(oc2[[1]])
  points(oc2[[1]], col=3)
}


