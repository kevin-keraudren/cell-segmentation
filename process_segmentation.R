#!/usr/bin/env Rscript

# load library + custom functions
library('EBImage', quiet=TRUE)
#source( file.path( Sys.getenv( 'R_LIBS' ), 'myfunctions.R' ) )
source( 'myfunctions.R' )

# read parameter file
parameter.file = commandArgs(trailingOnly=TRUE)[1]
source(parameter.file)

HECD1.images = Sys.glob(file.path(opt$folder,'HECD1',paste('*.',opt$format,sep='')))

if (length(HECD1.images)==0) {
print( file.path(opt$folder,'HECD1','*.',opt$format))
  print( '*** ERROR ***')
  print( 'No HECD1 images found' )
  print( 'Parameters are:')
  print( opt )
}

first.done = FALSE
for (hecd1.file in HECD1.images) {
  opt$n.debug = 1

  # get prefix for output file
  opt$tag = basename(hecd1.file)
  opt$tag = sub("\\.[^.]+$","",opt$tag, perl=TRUE)

  load(file.path(opt$outdir,paste(opt$tag,".RData",sep='')))
  data = list()

  # read raw data
  raw.hecd1 = readImage(hecd1.file)
  raw.actin = readImage(file.path(opt$folder,'ACTIN',basename(hecd1.file)))
  raw.dapi = readImage(file.path(opt$folder,'DAPI',basename(hecd1.file)))

  # size of images
  opt$size = dim(raw.hecd1)[1:2]
  
  # File id, cell id
  nb.cells = max(cells)
  data$file.id = rep(hecd1.file, nb.cells)
  data$cell.id = 1:nb.cells

  # Image border
  hf = hullFeatures(cells)
  data$touching.image.borders = as.integer(hf[,'g.edge'] > 0)

  # Background
  data$touching.background = rep(0, nb.cells)
  if (sum(background) != 0) {
    oc =ocontour(cells,external=TRUE)
    d = imageData(cells)
    for (i in 1:nb.cells ) {
      t = table( d[ oc[[i]] ] )
      if ("0" %in% names(t)) {
        data$touching.background[i] = 1
      }
    }
  }

  # Cell size and HECD1 in whole cell
  mvar = meanvar(cells,raw.hecd1)
  data$mean.hecd1.whole = mvar[,'mean']
  data$var.hecd1.whole = mvar[,'var']
  data$size.cell = mvar[,'size']

  # Actin in whole cell
  mvar = meanvar(cells,raw.actin)
  data$mean.actin.whole = mvar[,'mean']
  data$var.actin.whole = mvar[,'var']

  # Nuclei size
  mvar = meanvar(nuclei,raw.hecd1)
  data$size.nuclei = mvar[,'size']

  # Junctions
  cells.separation = separate(cells, kern = makeBrush(opt$analyse.erode, shape = "disc"))
  global.cytoplasm = cells.separation | (raw.hecd1 < opt$analyse.wall.threshold)

  ## Constant width junction
  walls.separate = cells
  walls.separate[cells.separation > 0] = 0
  cytoplasm.separate = cells
  cytoplasm.separate[cells.separation == 0] = 0
  cytoplasm.separate[nuclei > 0] = 0

  if (opt$debug) {
    img = saturate(raw.hecd1,opt$debug.saturate)
    segmentation = rgbImage(red=img, green=img, blue=img + saturate(raw.dapi,opt$debug.saturate))
    segmentation = paintObjects(cells.separation, segmentation, col=rgb(1,0,0)) 
    segmentation = paintObjects(walls.separate > 0, segmentation, col="cyan") 
    segmentation = paintObjects(nuclei, segmentation, col=rgb(0,1,0))
    segmentation = paintObjects(background, segmentation, col=rgb(1,0,1)) 
    segmentation = addNumbers( segmentation, cells)
    opt = debugImage( segmentation, 'constant_junction',opt)
  }
  
  mvar = meanvar(walls.separate,raw.hecd1)
  data$mean.hecd1.constant.junction = mvar[,'mean']
  data$var.hecd1.constant.junction = mvar[,'var']
  data$size.constant.junction = mvar[,'size']

  mvar = meanvar(cells.separation,raw.hecd1)
  data$mean.hecd1.constant.cytoplasm = mvar[,'mean']
  data$var.hecd1.constant.cytoplasm = mvar[,'var']
  data$size.constant.cytoplasm = mvar[,'size']

  ## Thresholded junction
  walls = cells
  walls[global.cytoplasm > 0] = 0
  cytoplasm = cells
  cytoplasm[global.cytoplasm == 0] = 0
  cytoplasm[nuclei > 0] = 0

  if (opt$debug) {
    segmentation = rgbImage(red=img, green=img, blue=img + saturate(raw.dapi,opt$debug.saturate))
    segmentation = paintObjects(cytoplasm, segmentation, col=rgb(1,0,0)) 
    segmentation = paintObjects(walls > 0, segmentation, col="cyan") 
    segmentation = paintObjects(nuclei, segmentation, col=rgb(0,1,0))
    segmentation = paintObjects(background, segmentation, col=rgb(1,0,1)) 
    segmentation = addNumbers( segmentation, cells)
    opt = debugImage( segmentation, 'thresholded_junction',opt)
  }
  
  mvar = meanvar(walls,raw.hecd1)
  data$mean.hecd1.thresholded.junction = mvar[,'mean']
  data$var.hecd1.thresholdedjunction = mvar[,'var']
  data$size.thresholded.junction = mvar[,'size']

  mvar = meanvar(cytoplasm,raw.hecd1)
  data$mean.hecd1.thresholded.cytoplasm = mvar[,'mean']
  data$var.hecd1.thresholded.cytoplasm = mvar[,'var']
  data$size.thresholded.cytoplasm = mvar[,'size']

  # Corner-to-corner
  contours = Image(0,dim=opt$size)
  oc = ocontour(cells,external=TRUE)
  for ( i in 1:nb.cells ) {
    contours[ oc[[i]] ] = contours[ oc[[i]] ] + 1
  }
  corners = contours > 1
  oc2 = list()
  length(oc2) = nb.cells
  oc = ocontour(cells,external=TRUE)
  hf = hullFeatures(cells)
  perimeter.ratio = c()
  for ( i in 1:nb.cells ) {
    if ( ! data$touching.image.borders[i]
        && sum(corners[oc[[i]]]) > 3 ) { # we exclude the round cells
      oc2[[i]] = oc[[i]][ corners[oc[[i]]],]
      oc2[[i]] = rbind(oc2[[i]],oc2[[i]][1,])
      perimeter.ratio[i] = perimeter(oc[[i]]) / perimeter( oc2[[i]] )
    }
    else {
      oc2[[i]] = array(c(0,0),dim=c(1,2))
      perimeter.ratio[i] = 1
    }
  }

  data$perimeter.ratio = perimeter.ratio

  if (opt$debug) {
    segmentation = rgbImage(red=img,green=img, blue=img+saturate(raw.dapi,opt$debug.saturate))
    segmentation = paintObjects(cells, segmentation, col=rgb(1,0,0)) # blue
    segmentation = paintObjects(nuclei, segmentation, col=rgb(0,1,0)) # red
    segmentation = paintObjects(background, segmentation, col=rgb(1,0,1)) # violet
    segmentation = drawPolyline(segmentation,oc2,stroke.color="midnightblue")
    segmentation = paintObjects(dilate(corners),segmentation,col=c("slateblue","slateblue") )
    segmentation = addNumbers( segmentation, cells)
    opt = debugImage( segmentation, 'corner-to-corner',opt)
  }
  
  # Write data into CSV file
  if (! first.done) {
    write.table(data,file=file.path(opt$outdir,"data.txt"),col.names=TRUE,row.names=FALSE,append=FALSE,
  sep='\t')
    first.done = TRUE
  }
  else {
    write.table(data,file=file.path(opt$outdir,"data.txt"),col.names=FALSE,row.names=FALSE,append=TRUE,
  sep='\t')
  }

  if (opt$debug) {
    warnings()
  }

  if (opt$only.one) {
    quit()
  }
  
}
