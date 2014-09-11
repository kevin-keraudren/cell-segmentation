\name{magickFilter}

\alias{maxFilter}
\alias{minFilter}
\alias{meanFilter}
\alias{medianFilter}
\alias{gradientFilter}
\alias{nonpeakFilter}
\alias{stdFilter}
\alias{modeFilter}


\title{ImageMagick StatisticImage filters}

\description{
  ImageMagick StatisticImage filters: min, max, mean, median, gradient,
  nonpeak, standard deviation, mode.
}

\usage{
  minFilter(x,w=3,h=3)
  maxFilter(x,w=3,h=3)
  meanFilter(x,w=3,h=3)
  medianFilter(x,w=3,h=3)
  gradientFilter(x,w=3,h=3)
  nonpeakFilter(x,w=3,h=3)
  stdFilter(x,w=3,h=3)
  modeFilter(x,w=3,h=3)

}

\arguments{
  \item{x}{An \code{Image} object or an array.}
  
  \item{w}{A numeric value for the width of the pixel neighbourhood.}

  \item{h}{A numeric value for the height of the pixel neighbourhood.}
  
}

\value{
  An \code{Image} object or an array, containing the filtered version
  of \code{x}.
}

\details{
  These filters only use a box window.
}

\seealso{
  \code{filter2}
}

\references{
    \emph{ImageMagick StatisticImage}: \url{http://www.imagemagick.org/api/effect.php#StatisticImage}.
}

\author{
  Kevin Keraudren, \email{kevin.keraudren10@imperial.ac.uk}
}

\examples{
   x = readImage(system.file("images", "lena.gif", package="EBImage"))
   if (interactive()) display(x)

   y = minFilter(x)
   if (interactive()) display(y, title='minFilter')

   y = maxFilter(x)
   if (interactive()) display(y, title='maxFilter')

   y = meanFilter(x)
   if (interactive()) display(y, title='meanFilter')

   y = medianFilter(x)
   if (interactive()) display(y, title='medianFilter')

   y = gradientFilter(x)
   if (interactive()) display(y, title='gradientFilter')

   y = nonpeakFilter(x)
   if (interactive()) display(y, title='nonpeakFilter')

   y = stdFilter(x)
   if (interactive()) display(y, title='stdFilter')

   y = modeFilter(x)
   if (interactive()) display(y, title='modeFilter')

}

