\name{separate}

\alias{separate}

\title{Separate touching objects}

\description{
  Functions to perform morphological erosion to separate touching objects.
}

\usage{
  separate(x, kern = makeBrush(5, shape = "diamond")) 
}

\arguments{
  \item{x}{An \code{Image} object or an array containing labelled regions.}
  
  \item{kern}{An \code{Image} object or an array, containing the
    structuring element. \code{kern} is considered as a binary image, whose
    pixels of value 0 are considered as background ones and other pixels as
    foreground ones.}
}

\value{
  Returns the
  transformed \code{Image} object or array,  after separating the
  touching objects.
}

\author{
  Kevin Keraudren, \email{kevin.keraudren10@imperial.ac.uk}
}

\examples{
  x = readImage(system.file('images', 'shapes.png', package='EBImage'))
  x = x[110:512,1:130]
  if (interactive()) display(x, title='Binary')
  y = distmap(x)
  if (interactive()) display(normalize(y), title='Distance map')
  w = watershed(y)
  if (interactive()) display(normalize(w), title='Watershed')	
  z = separate(w)
  if (interactive()) display(normalize(z), title='Separate')
}


