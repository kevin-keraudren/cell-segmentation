#!/usr/bin/env Rscript

# load library + custom functions
library('EBImage', quiet=TRUE)
#source( file.path( Sys.getenv( 'R_LIBS' ), 'myfunctions.R' ) )
source( 'myfunctions.R' )

# read parameter file
parameter.file = commandArgs(trailingOnly=TRUE)[1]
source(parameter.file)

opt$first.watergrow.kernel = makeBrush(opt$first.watergrow.size, shape='disc')
opt$first.watergrow.min.votes = as.integer(sum(opt$first.watergrow.kernel)/opt$first.watergrow.ratio)
opt$second.watergrow.kernel = makeBrush(opt$second.watergrow.size, shape='disc')
opt$second.watergrow.min.votes = as.integer(sum(opt$second.watergrow.kernel)/opt$second.watergrow.ratio)

# create ouput_dir and temp_dir if they do not exist
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
dir.create(opt$tmp.dir, showWarnings=FALSE, recursive=TRUE)


HECD1.images = Sys.glob(file.path(opt$folder,'HECD1',paste('*.',opt$format,sep='')))

if (length(HECD1.images)==0) {
print( file.path(opt$folder,'HECD1','*.',opt$format))
  print( '*** ERROR ***')
  print( 'No HECD1 images found' )
  print( 'Parameters are:')
  print( opt )
}


for (hecd1.file in HECD1.images) {
  opt$n.debug = 1

  # get prefix for output file
  opt$tag = basename(hecd1.file)
  opt$tag = sub("\\.[^.]+$","",opt$tag, perl=TRUE)

  # read raw data
  raw.hecd1 = readImage(hecd1.file)
  raw.actin = readImage(file.path(opt$folder,'ACTIN',basename(hecd1.file)))
  raw.dapi = readImage(file.path(opt$folder,'DAPI',basename(hecd1.file)))

  # size of images
  opt$size = dim(raw.hecd1)[1:2]

  # Preprocessed data
  if (opt$do.preprocessing) {
  clean.hecd1 = ced(raw.hecd1,opt)
  clean.actin = eed(raw.actin,opt)
  clean.dapi = eed(raw.dapi,opt)
  }
  else {
  clean.hecd1 = raw.hecd1
  clean.actin = raw.actin
  clean.dapi = raw.dapi  
  }

  # compute Sobel edge map and perform morphological closing
  edge.map = sobel(saturate(clean.hecd1,opt$sobel.saturate))
  if (opt$debug) opt = debugImage( edge.map,'Sobel', opt)
  edge.map = closing(edge.map, kern=makeBrush(opt$edge.map.kernel.size, shape='disc'), binary=FALSE)
  if (opt$debug) opt = debugImage( edge.map, 'Edge-map',opt)

  # threshold nuclei
  thresholded.dapi = thresh(clean.dapi, opt$dapi.adaptive.thresholding.size,opt$dapi.adaptive.thresholding.size, opt$dapi.adaptive.thresholding.offset)
  if (opt$debug) opt = debugImage( thresholded.dapi, 'thresholded_dapi',opt)
  thresholded.dapi = fillHull( thresholded.dapi )
  if (opt$debug) opt = debugImage( thresholded.dapi, 'thresholded_dapi',opt)
  if ( opt$split.nuclei ) {
	  distmap.dapi = distmap( thresholded.dapi )
	  watershed.dapi = watershed(distmap.dapi)
	  if (opt$debug) opt = debugImage( watershed.dapi,'watershed_dapi',opt)
  } else {
	  watershed.dapi = bwlabel( thresholded.dapi)
  }

  # threshold background
  intensity = pmax( as.vector(clean.hecd1), as.vector(clean.actin),
    as.vector(clean.dapi))
  intensity = Image(intensity, opt$size )
  background = intensity < opt$background.threshold
  background = bwlabel( background )
  if (opt$debug) {
	  segmentation = rgbImage(red=saturate(raw.hecd1,opt$debug.saturate),
		  green=saturate(raw.actin,opt$debug.saturate), blue=saturate(raw.dapi,opt$debug.saturate))
	  segmentation = paintObjects(background, segmentation, col=rgb(1,0,1)) # violet
	opt = debugImage( segmentation, 'initial_background',opt)
  }


  # remove small objects from background and from not.background
  hf = hullFeatures( background )    
  background = rmObjects(background, which(hf[,"g.s"] < opt$small.objects ))    
  not.background = bwlabel(background == 0)
  hf = hullFeatures( not.background )
  not.background = rmObjects(not.background, which(hf[,"g.s"] < opt$small.objects ))
  background = not.background == 0

  # dilate background to avoid leaks
  if (opt$background.dilate > 0) {
    background = dilate(background, kern=makeBrush(opt$background.dilate,'disc'))
    background = background > 0
  }
  if (opt$debug) {
	  segmentation = rgbImage(red=saturate(raw.hecd1,opt$debug.saturate),
		  green=saturate(raw.actin,opt$debug.saturate), blue=saturate(raw.dapi,opt$debug.saturate))
	  segmentation = paintObjects(background, segmentation, col=rgb(1,0,1)) # violet
	opt = debugImage( segmentation, 'initial_background',opt)
  }

  # remove small nuclei
  hf = hullFeatures(watershed.dapi)
  watershed.dapi = rmObjects(watershed.dapi, which(hf[,"g.s"] < opt$small.objects ))
  watershed.dapi = reenumerate(watershed.dapi)

  # get geometric center of nuclei
  hf = hullFeatures(watershed.dapi)
  nuclei.centers = Image(dim=opt$size)
  nuclei.centers[hf[,c("g.x","g.y")]] = 1
  tmp.nuclei.centers = dilate(nuclei.centers,kern=opt$first.watergrow.kernel)
  nuclei.centers = watershed.dapi
  nuclei.centers[tmp.nuclei.centers == 0] = 0

  # ignore background in edge map
  tmp.edge.map = edge.map
  tmp.edge.map[background] = -2

  # first watergrow using edge map
  watergrow.result = watergrow(tmp.edge.map, nuclei.centers,
    min.water=0,kern=opt$first.watergrow.kernel, min.votes=opt$first.watergrow.min.votes)

  watergrow.result = fillHull(watergrow.result)
  background = watergrow.result == 0

  # remove wrongly classified background
  if (sum(background) != 0) {
    bg = reenumerate(bwlabel(background),start=max(watergrow.result)+1)
    bgwatergrow = overlay(watergrow.result,bg)
    oc =ocontour(bgwatergrow,external=TRUE)
    data = imageData(bgwatergrow)
    for (i in ((max(watergrow.result)+1):max(bg)) ) {
      t = table( data[ oc[[i]] ] )
      t[[ as.character(i) ]] = 0
      if (max(t) / sum(t) > 0.7) {
        watergrow.result[bg == i] = as.integer(names(which.max(t)))
      }
    }
    background =  watergrow.result == 0
  }

  if (opt$debug) {
	segmentation = rgbImage(red=clean.hecd1,
		green=clean.hecd1, blue=clean.hecd1 + clean.dapi)
	segmentation = paintObjects(nuclei.centers, segmentation, col=rgb(0,1,0)) 
	segmentation = paintObjects(watergrow.result, segmentation, col=rgb(1,0,0))
	segmentation = paintObjects(background, segmentation, col=rgb(1,0,1)) 
	opt = debugImage( segmentation, 'first_watershed',opt)
  }

  if ( opt$refinement.step ) {

	# separate cells before the last watergrow (we shrink the cells)
	watergrow.result = separate(watergrow.result,kern=makeBrush(opt$shrink.cells, shape='disc'))
	watergrow.result = keep.bigger.part(watergrow.result)

	# remove cells that are too small (likely multi-nuclei cells)
	hf = hullFeatures(watergrow.result)
	watergrow.result = rmObjects(watergrow.result, which(hf[,"g.s"] < opt$small.objects ))
	watergrow.result = reenumerate(watergrow.result)
	watergrow.result = fillHull(watergrow.result)

	# run a watergrow, this time on the raw HECD1
	hecd1 = raw.hecd1
	hecd1[background != 0] = -2
	watergrow.result = watergrow(hecd1, watergrow.result, min.water=0,kern=opt$second.waergrow.kernel, min.votes=opt$second.watergrow.min.votes )

	# we will remove cells:
	# - that are too small
	# - when nuclei are not separated by "walls"
	cells.separation = separate(watergrow.result, kern = makeBrush(5, shape = "disc"))
	mvar = meanvar(watergrow.result, raw.hecd1)
	values = mvar[, "mean"] + 4.0 * mvar[, "var"]
	cell.thresh = imageReplace(watergrow.result, values, do.zero = FALSE)
	global.cytoplasm = cells.separation | (raw.hecd1 < cell.thresh)
	walls = watergrow.result > 0 & global.cytoplasm == 0

	mask = Image(dim=opt$size,colormode='grayscale')
	mask = paintObjects(watergrow.result,mask,col='white')
	candidates = mask & watershed.dapi & ! walls
	candidates = closing(candidates,kern = makeBrush(11, shape = "disc"))
	candidates = bwlabel(candidates)
	hf = hullFeatures(candidates)
	candidates = rmObjects(candidates,which(hf[,'g.s'] < opt$min.nuclei.size.for.merging))

	if (max(candidates) > 0 ) {
	candidates = dilate(candidates,kern = makeBrush(3, shape = "box"))
	candidates = bwlabel(candidates)
	data = imageData(watergrow.result)
	for (i in (1:max(candidates)) ) {
	  f = factor( data[ candidates==i ] )
	  l = levels(f)
	  l1 = as.numeric(l[[1]])
	  for (n in levels(f)) {
		
	   n.i = as.numeric(n)
	   if (n.i == 0) next # background
	   if (n.i == l1) next
	   watergrow.result[watergrow.result == n.i] = l1
	  }
	}
	watergrow.result = reenumerate(watergrow.result)
	}

	# do it again for cells
	cells.separation = separate(watergrow.result, kern = makeBrush(5, shape = "disc"))
	mvar = meanvar(watergrow.result, raw.hecd1)
	values = mvar[, "mean"] + 4.0 * mvar[, "var"]
	cell.thresh = imageReplace(watergrow.result, values, do.zero = FALSE)
	global.cytoplasm = cells.separation | (raw.hecd1 < cell.thresh)
	walls = watergrow.result > 0 & global.cytoplasm == 0

	mask = Image(dim=opt$size,colormode='grayscale')
	mask = paintObjects(watergrow.result,mask,col='white')
	candidates = mask & ! walls
	candidates = closing(candidates,kern = makeBrush(5, shape = "disc"))
	candidates = bwlabel(candidates)
	hf = hullFeatures(candidates)
	candidates = rmObjects(candidates,which(hf[,'g.s'] < opt$min.cell.size.for.merging))

	if (max(candidates) > 0 ) {
	candidates = dilate(candidates,kern = makeBrush(3, shape = "box"))
	candidates = bwlabel(candidates)
	data = imageData(watergrow.result)
	for (i in (1:max(candidates)) ) {
	  f = factor( data[ candidates==i ] )
	  l = levels(f)
	  l1 = as.numeric(l[[1]])
	  for (n in levels(f)) {
		
	   n.i = as.numeric(n)
	   if (n.i == 0) next # background
	   if (n.i == l1) next
	   watergrow.result[watergrow.result == n.i] = l1
	 }
	}
	watergrow.result = reenumerate(watergrow.result)
	}
	
    if (opt$debug) {
	  segmentation = rgbImage(red=clean.hecd1,
	  	green=clean.hecd1, blue=clean.hecd1 + clean.dapi)
	  segmentation = paintObjects(nuclei.centers, segmentation, col=rgb(0,1,0)) 
	  segmentation = paintObjects(watergrow.result, segmentation, col=rgb(1,0,0)) 
	  segmentation = paintObjects(background, segmentation, col=rgb(1,0,1)) 
	  opt = debugImage( segmentation, 'refinement_step',opt)
    }	

  }

  # separate cells before the last watergrow (we shrink the cells)
  watergrow.result = separate(watergrow.result,kern=makeBrush(opt$shrink.cells, shape='disc'))
  watergrow.result = keep.bigger.part(watergrow.result)

  # remove cells that are too small (likely multi-nuclei cells)
  hf = hullFeatures(watergrow.result)
  watergrow.result = rmObjects(watergrow.result, which(hf[,"g.s"] < opt$small.objects ))
  watergrow.result = reenumerate(watergrow.result)
  watergrow.result = fillHull(watergrow.result)

  # run a watergrow, this time on the raw HECD1
  hecd1 = raw.hecd1
  hecd1[background != 0] = -2
  watergrow.result = watergrow(hecd1, watergrow.result, min.water=0,kern=opt$second.watergrow.kernel, min.votes=opt$second.watergrow.min.votes )

  # redefine the nuclei
  cells = watergrow.result
  nuclei = cells
  nuclei[watershed.dapi == 0] = 0
  nuclei = separate(nuclei,kern=makeBrush(3, shape = "diamond"))

  # output image of segmentation
  img = saturate(raw.hecd1,opt$debug.saturate)
  segmentation = rgbImage(red=img, green=img, blue=img + saturate(raw.dapi,opt$debug.saturate))
  segmentation = paintObjects(nuclei, segmentation, col=rgb(0,1,0))
  segmentation = paintObjects(cells, segmentation, col=rgb(1,0,0))
  segmentation = paintObjects(background, segmentation, col=rgb(1,0,1))

  hf = hullFeatures(cells)
  xy = hf[, c('g.x', 'g.y')]
  labels = as.character(1:max(cells))
  font = drawfont(weight=600, size=20)
  segmentation = drawtext(segmentation, xy=xy, labels=labels , font=font, col="white")

  writeImage( segmentation, file.path(opt$outdir,basename(hecd1.file)))


  # save segmentation
  save(cells,nuclei,background,file=paste(opt$outdir,'/',opt$tag,'.RData',sep=''))

  if (opt$debug) {
    warnings()
  }

  if (opt$only.one) {
    quit()
  }
}


