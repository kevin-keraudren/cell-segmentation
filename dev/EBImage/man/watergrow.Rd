\name{watergrow}

\alias{watergrow}

\concept{object detection}

\title{Seeded watershed}

\description{
Find boundaries between adjacent regions in an image, where seeds
have been already identified in the individual regions to be
segmented.
}

\usage{
  watergrow( map, labels,
  min.water=min(map), max.water=max(map),
  step=(quantile(map,0.98) - quantile(map,0.2))/20,
  kern=makeBrush(5, shape='diamond'),
  min.votes=4
  )
}

\arguments{
   \item{map}{An \code{Image} object or an array, containing the image to segment.}

   \item{labels}{An \code{Image} object or an array, containing the seeding objects of
   the already identified regions.}	

   \item{min.water}{Starting water level, pixels below this level will
   be ignored.}

   \item{max.water}{Maximum water level, pixels above this level will be
   ignored.}

   \item{step}{Amount of water added at each step.}

  \item{kern}{Kernel used to count the number of votes for a pixel to be
  assigned a label.}

  \item{min.votes}{Minimum number of votes required for a pixel to be
    assigned a label.}
}

\value{
  An \code{Image} object or an array, containing the labelled objects.
}

\details{
  The water starts at \code{min.level} and raises until it reaches
  \code{max.level}, adding \code{step} water at each step. 
}

\seealso{ \code{\link{bwlabel}}, \code{\link{watershed}}, \code{\link{propagate}}
}

\examples{

## Watergrow avoiding leaks
x = readImage(system.file('images', 'watergrow-test.png', package='EBImage'))
nuclei = x[,,1]
colorMode(nuclei) = 'Grayscale'
nuclei = bwlabel(nuclei)
borders = x[,,2]
colorMode(borders) = 'Grayscale'
segmentation = watergrow(borders, nuclei, kern = makeBrush(11, shape
 = "disc"), min.votes = 40,max.water=0.9)
if (interactive()) display(normalize(segmentation), title='watergrow segmentation')

}

\author{
  Kevin Keraudren, \email{kevin.keraudren10@imperial.ac.uk}
}



