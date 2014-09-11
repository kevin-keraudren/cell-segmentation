\name{drawPolyline}

\alias{drawPolyline}

\title{Draw a polyline on an image.}

\description{
  Draw a polyline on an image (serie of segments).
}

\usage{
  drawPolyline(img, xy, stroke.color = "white", stroke.width = 1, 
    fill.color = "none", antialias = FALSE)
}

\arguments{
\item{img}{An \code{Image} object or an array.}  
\item{xy}{list of vectors of coordinates}
  
\item{stroke.color}{A numeric or a character string specifying the color
  of the polyline.}

\item{stroke.width}{A numeric specifying the width
    of the polyline.}
  
\item{fill.color}{A numeric or a character string specifying the color
  of to fill the area defined by the polyline. \code{"none"} stands for
  no filling.}
  
\item{antialias}{A logical indicating whether to turn antialiasing on.}
}

\value{
  An \code{Image} object or an array, containing the transformed version
  of \code{img}.
}

\author{
  Kevin Keraudren, 2011
}

\examples{

img = readImage(system.file("images", "shapes.png", package="EBImage"))
img = img[110:512,1:130]
labels = bwlabel(img)
oc = ocontour(labels)
poly = list()
for ( i in 1:max(labels) ) {
  x = filter(oc[[i]][,1],rep(1,21)/21,method='convolution',circular=TRUE)
  y = filter(oc[[i]][,2],rep(1,21)/21,method='convolution',circular=TRUE)
  p = cbind(x,y)
  p = rbind(p,p[1,])
  poly[[i]] = p
}
output = rgbImage(red=img,green=img,blue=img)
output = drawPolyline(output,poly,stroke.color='red')
if (interactive()) display(output)

}


