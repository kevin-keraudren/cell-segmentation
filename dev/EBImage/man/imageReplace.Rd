\name{imageReplace}

\alias{imageReplace}

\title{Replace labels by values}

\description{
  Replaces labels in a mask by values from a list.
}

\usage{
  imageReplace(x, colors, do.zero = FALSE) 
}

\arguments{
  \item{x}{An \code{Image} object or an array containing object masks.
    Object masks are sets of pixels with the same unique integer value.}
    
  \item{colors}{List of values that will replace the labels.}
  
  \item{do.zero}{Logical value indicating whether the first value in the
    \code{colors} list stands for the label 0 or the label 1.}
}

\value{
  An \code{Image} object or an array with the labels replaced by the
  corresponding values.
}


\seealso{
  \code{\link{colorize}}
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

  ## moments
  m = moments(mask, x)

  img = imageReplace(mask[,,1],m[[1]][,'m.ecc'])
  if (interactive()) display(colorize(img),
      title='Nuclei colored by their eccentricity (red is high)')
}


