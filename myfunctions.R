sobel = function( img ) {
  data = c(1,0,-1,2,0,-2,1,0,-1)
  kernelX = matrix(data, 3, 3, byrow = TRUE)
  kernelY = matrix(data, 3, 3, byrow = FALSE)
  x = imageConvolve(img, kernelX)
  y = imageConvolve(img, kernelY)
  return( sqrt(x^2 + y^2) )
}

gradient = function( img ) {
  data = c(-1,0,-1,1,0,-1,1,0,-1)
  kernelX = matrix(data, 3, 3, byrow = TRUE)
  kernelY = matrix(data, 3, 3, byrow = FALSE)
  x = magickConvolve(img, kernelX)
  y = magickConvolve(img, kernelY)
  return( sqrt(x^2 + y^2) )
}

keep.bigger.part = function(mask) {
  buffer = bwlabel(mask) # indexing the children
  hf = hullFeatures(buffer)
  cm = cmoments(buffer,mask)
  parents = cm[,'m.int']/cm[,'m.pxs'] # index of the parents
  indexes = c() # indices to delete
  for (i in 1:max(mask)) {
    index =  which(hf[,"g.s" ] == max(hf[ parents == i,"g.s" ]) & parents == i)
    indexes = c(indexes,index)
  }
  buffer = rmObjects(buffer, indexes)
  mask[ buffer != 0 ] = 0
  return(mask)
}

saturate = function(img,q) {
  t = quantile(img,q)
  img[img > t] = t
  return(normalize(img))
}

saturate2 = function(img,q1,q2) {
  t = quantile(img,q1)
  img[img < t] = t
  t = quantile(img,q2)
  img[img > t] = t  
  return(normalize(img))
}

overlay = function(below, above) {
  below[above != 0] = 0
  return(below + above)
}

remove.border.objects = function(mask) {
  hf = hullFeatures(mask)
  return( rmObjects( mask, which( hf[,'g.edge'] > 0 ) ) )
}

eed = function(img,opt) {
  cparam = opt$eed.cparam
  iteration = opt$eed.iteration
  tmp.name = tempfile(pattern = "eed", fileext = "", tmpdir=opt$tmp.dir)
  file.create( tmp.name )
  tmp.1.png = paste(tmp.name, '.png',sep='')
  tmp.1.pgm = paste(tmp.name, '.pgm',sep='')
  writeImage(img, tmp.1.png,quality=100)
  if (opt$OS == 'Linux') convert.prog = 'convert'
  else convert.prog = 'setenv.bat && convert'
  system(paste(convert.prog, ' ', tmp.1.png, ' ', tmp.1.pgm, sep=''))
  if (opt$OS == 'Linux') eed.prog = 'eed'
  if (opt$OS == 'Win32') eed.prog = 'setenv.bat && eed_win32b.exe'
  if (opt$OS == 'Win64') eed.prog = 'setenv.bat && eed_win64b.exe'
  system(paste(eed.prog, ' ',tmp.1.pgm, ' ',cparam,' ',iteration,' ',iteration,' ',tmp.name,sep=''))
  tmp.2.pgm = paste( tmp.name, '_', sprintf('%04d',iteration), '.pgm',sep='')
  out = readImage( tmp.2.pgm )
  file.remove( tmp.1.pgm )
  file.remove( tmp.1.png )
  if (! opt$debug) file.remove( tmp.2.pgm )
  file.remove( tmp.name )
  return( out )
}

ced = function(img,opt) {
  rho = opt$ced.rho
  iteration = opt$ced.iteration
  tmp.name = tempfile(pattern = "ced", fileext = "", tmpdir=opt$tmp.dir)
  file.create( tmp.name )
  tmp.1.png = paste(tmp.name, '.png',sep='')
  tmp.1.pgm = paste(tmp.name, '.pgm',sep='')
  writeImage(img, tmp.1.png,quality=100)
  if (opt$OS == 'Linux') convert.prog = 'convert'
  else convert.prog = 'setenv.bat && convert'
  system(paste(convert.prog, ' ', tmp.1.png, ' ', tmp.1.pgm, sep=''))
  out = sprintf('%04d',iteration)
  if (opt$OS == 'Linux') ced.prog = 'ced'
  if (opt$OS == 'Win32') ced.prog = 'setenv.bat && ced_win32b.exe'
  if (opt$OS == 'Win64') ced.prog = 'setenv.bat && ced_win64b.exe'
  system(paste(ced.prog, ' ',tmp.1.pgm, ' ',rho,' ',iteration,' ',iteration,' ',tmp.name,sep=''))
  tmp.2.pgm = paste( tmp.name, '_', sprintf('%04d',iteration), '.pgm',sep='')
  out = readImage( tmp.2.pgm )
  file.remove( tmp.1.pgm )
  file.remove( tmp.1.png )
  if (! opt$debug) file.remove( tmp.2.pgm )
  file.remove( tmp.name )
  return( out )
}

debugImage = function(img,name, opt) {
  if (colorMode(img) == 'Grayscale') img = normalize(img)
  writeImage( img, file.path(opt$tmp.dir, paste('debug_',opt$tag,'_',opt$n.debug, '_',name,'.png',sep='')))
  opt$n.debug = opt$n.debug + 1
  return(opt)
}

addNumbers = function( img, objects ) {
  hf = hullFeatures(objects)
  xy = hf[, c('g.x', 'g.y')]
  labels = as.character(1:max(objects))
  font = drawfont(weight=600, size=20)
  img = drawtext(img, xy=xy, labels=labels , font=font, col="white")
  return(img)
 }
  
perimeter = function(poly) {
  i <- 1:(length(poly[,1]) - 1)
  return( sum(sqrt( (poly[i,1] -poly[i+1,1])^2 + (poly[i,2] - poly[i + 1,2])^2) ) )
}


