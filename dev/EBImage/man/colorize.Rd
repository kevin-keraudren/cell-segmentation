\name{colorize}

\alias{colorize}

\title{Colormap}

\description{
  Apply a colormap to an image.
}

\usage{
  colorize(img, colormap = "jet") 
}

\arguments{
  \item{img}{A grayscale \code{Image} object or an array.}
  
  \item{colormap}{A character value for the colormap. Currently accepted
  values are "jet" or "rainbow".}
  
}

\value{
  A color \code{Image} object or an array.
}

\author{
  Kevin Keraudren, \email{kevin.keraudren10@imperial.ac.uk}
}

\examples{
   x = readImage(system.file("images", "lena.gif", package="EBImage"))
   if (interactive()) display(colorize(x))
}

