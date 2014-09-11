\name{imageConvolve}

\alias{imageConvolve}

\title{2D Convolution Filter}

\description{
  Convolves an image with the given kernel.
}

\usage{
  imageConvolve(x, filter)
}

\arguments{
  \item{x}{An \code{Image} object or an array.}

  \item{filter}{An \code{Image} object or an array, with odd spatial
    dimensions. Must contain only one frame.}
}

\value{
  An \code{Image} object or an array, containing the filtered version
  of \code{x}.
}

\details{
  Linear filtering is useful to perform low-pass filtering (to blur
  images, remove noise...) and high-pass filtering (to detect
  edges, sharpen images). The function \code{makeBrush} is useful to
  generate filters.

  Data is set to 0 around borders.
  
  If \code{x} contains multiple franes, the filter will be applied one each frame.
  
}

\seealso{
  \code{\link{makeBrush}}, \code{\link{convolve}}, \code{\link{blur}},\code{\link{minFilter}}, \code{\link{filter2}}
}

\author{Kevin Keraudren, \email{kevin.keraudren10@imperial.ac.uk}}

\examples{
  x = readImage(system.file("images", "lena-color.png", package="EBImage"))
  if (interactive()) display(x, title='Lena')

  # Sobel filter for edge detection
  sobel = function( img ) {
    data = c(1,0,-1,2,0,-2,1,0,-1)
    kernelX = matrix(data, 3, 3, byrow = TRUE)
    kernelY = matrix(data, 3, 3, byrow = FALSE)
    x = imageConvolve(img, kernelX)
    y = imageConvolve(img, kernelY)
    return( sqrt(x^2 + y^2) )
  }

  y = sobel(x)
  if (interactive()) display(normalize(y), title='Sobel lena')

}
