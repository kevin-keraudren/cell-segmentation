#######################
### Parameter files ###
#######################
#
# lines beginning by a '#' are comments
# values TRUE or FALSE do not take quotes
#
# Written by Kevin Keraudren, 04/11/2011
# kevin.keraudren10@imperial.ac.uk
#

### GENERAL ###
opt = list() # do not change this line

opt$OS = 'Linux'  # either 'Win32', 'Win64' or 'Linux'

opt$folder = 'input_images' # This folder contains 3 subfolders: HECD1, DAPI and ACTIN
                              # In each subfolder, corresponding images need to
                              # have the same name 

opt$outdir = 'output_dir' # directory for resulting segmentation (*.png) and
                            # data.csv file (to be used with Excel)

opt$tmp.dir = 'tmp' # directory for temporary and debug output

opt$format = 'tif' # Filename extension of image files (do not include '.')

opt$debug = TRUE # for images of each step (useful to tweak parameters)
opt$only.one = FALSE # process only one image (useful to tweak parameters)
opt$debug.saturate = 1.0 # saturate channels in debugging output (1.0 does no
                         # saturation, 0.98 is good for dim images)

### Preprocessing step ###
opt$do.preprocessing = TRUE # skipping preprocessing allows to get a quick draft

# CED is an edge enhancing algorithm (used for HECD1 image)
# if in doubt, try rho=10, iteration=30
opt$ced.rho = 10
opt$ced.iteration = 30

# EED is a denoising algorithm which preserves edges
# if in doubt, try cparam=10, iteration=30
opt$eed.cparam = 10
opt$eed.iteration = 30


### Main processing (two-step watershed) ###
opt$edge.map.kernel.size = 21 # morphological closing for the edge map

opt$background.threshold = 0.08 # threshold for the background; thresholding is
                                # done on the maximum intensity of the 3 channels
opt$background.dilate = 21 # dilatation of the background to avoid watershed leaks

opt$dapi.adaptive.thresholding.size = 50 # size of the window for the adaptive
                                        # thresholding on the DAPI
opt$dapi.adaptive.thresholding.offset=0.05 # offset for the adaptive
                                        # thresholding on the DAPI

opt$small.objects = 800 # nuclei and cells below this size will be ignored

opt$first.watergrow.size = 31 # kernel size (must be odd)
opt$first.watergrow.ratio = 2.8 # ratio of the kernel which must already be in
                                # the region for a new pixel to be
                                # included; this is used to avoid watershed
                                # leaks

opt$shrink.cells = 9 # eroding cells before the second watergrow

opt$second.watergrow.size = 7 # this value should be odd and smaller than
                              # opt$first.watergrow.size 
opt$second.watergrow.ratio = 5 # this value should be greater than
                               # opt$first.watergrow.ratio 

# Advanced parameters
opt$sobel.saturate=1.0 # saturate the edge-map, useful in case of low contrast
                       # leave 1.0 for no saturation
opt$split.nuclei=FALSE # split adjacent nuclei
opt$refinement.step=FALSE # merge over-segmented cells (defines multi-nuclei cells)
opt$min.nuclei.size.for.merging = 20
opt$min.cell.size.for.merging = 80

### Analysis ###
opt$analyse.erode = 3 # erosion to define the junctions
opt$analyse.wall.threshold = 0.6 # threshold to define the junctions
