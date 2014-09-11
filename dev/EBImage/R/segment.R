## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
watershed = function (x, tolerance=1, ext=1) {
  validImage(x)
  if (colorMode(x)==TrueColor) stop("'x' must be an Image not in \'TrueColor\' color mode")
  tolerance = as.numeric(tolerance)
  if (tolerance < 0) stop( "'tolerance' must be non-negative" )
  ext = as.integer(ext)
  if (ext<1) stop( "'ext' must be a positive integer" )
  .Call("watershed", castImage(x), tolerance, ext, PACKAGE='EBImage')
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
watergrow = function ( map, labels,
  min.water=min(map), max.water=max(map),
  step=(quantile(map,0.98) - quantile(map,0.2))/20,
  kern=makeBrush(5, shape='diamond'),
  min.votes=4
  ) {
  validImage(map)
  if (colorMode(map)==TrueColor) stop("'map' must be an Image not in \'TrueColor\' color mode")
  min.water = as.numeric(min.water)
  max.water = as.numeric(max.water)
  step = as.numeric(step)
  min.votes = as.integer(min.votes)
  if (step==0) stop( "'step_size' must not be zero" )
  .Call("watergrow", castImage(map),castImage(labels),
        min.water, max.water,step,
        kern, min.votes,
        PACKAGE='EBImage')
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
propagate = function (x, seeds, mask=NULL, lambda=1e-4, ext, seed.centers) {
  validImage(x)
  checkCompatibleImages(x, seeds)
  
  if (!is.null(mask)) {
    checkCompatibleImages(x, mask)
    mask = castImage(mask)
  }

  lambda = as.numeric(lambda)
  if (lambda < 0.0) stop("'lambda' must be positive" )

  if (!missing(ext)) warning("'ext' is deprecated.")
  
  if (!missing(seed.centers)) warning("'seed.centers' is deprecated.")
  else seed.centers = FALSE
  
  if (seed.centers) {
    cm = hullFeatures(seeds)
    dimx = dim(seeds)
    nz = getNumberOfFrames(seeds, 'total')
    if (nz==1) cm = list(cm)
    index = lapply(cm, function(xy) floor(xy[,'g.x']) + floor(xy[,'g.y'])*dimx[1])
    for (i in 1:nz) index[[i]] = index[[i]] + (i-1)*dimx[1]*dimx[2]
    index = unlist(index)
    s = Image(0, dim=dim(seeds))
    s[index] = seeds[index]
    seeds = s
  }
  
  return(.Call( "propagate", castImage(x), castImage(seeds), mask, lambda, PACKAGE='EBImage'))
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ocontour = function(x, external=FALSE) {
  validImage(x)
  storage.mode(x)='integer'
  y = .Call('ocontour', x, as.integer(external), PACKAGE='EBImage')[-1]
  y = lapply(y, function(z) matrix(z, nc=2, byrow=TRUE))
  y
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bwlabel = function(x) {
  validImage(x)
  .Call("bwlabel", castImage(x), PACKAGE='EBImage')
}
