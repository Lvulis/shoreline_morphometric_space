### Script to convert occurrence into watermask, run OAM, and then extract shoreline shapefile.

##### Load libraries #####

library(ROAM)
library(stars)
library(raster)
library(sf)

# this script requires  `imager`, which is not loaded not loaded explicitly but is used for all image processing computations.
# `EBImage` is used for image visualization and some functions are taken from it,
# It was not used for image processing due to it not supporting 8-neighbor connectivity.#

#### Declare some constants and set up workspace ###
# Change to working directory

dirname <- "~/deltas/shoreline_morphometric_space/analysis/"

setwd(dirname)

# load utility functions
source("./scripts/utils.R")
dname = "Selenga"
pixres = 30 # Project to a 30 m spatial resolution image
setwd(paste0(dirname, dname))

#### Begin analysis ####
# Need to project GSW from WGS84 (EPSG:4326) to local coordinate system.
# Use local UTM zone and use a point near delta apex to estimate UTM zone
# Can skip this if revisiting delta and want to use different watermask threshold
# But don't need to rerun occurrence
write_occ = T
fn_oc_UTM <- "./occurrence_UTM.tiff"
if(write_occ) {
  # Load occurrence
  occurrence = read_stars("occurrence.tif")
  
  reference_dot <- st_read(paste0("./outlines/",dname,"_reference_dot_unprojected.shp"))
  
  lonlat <- st_coordinates(reference_dot) # Get lonlat
  inlet_node_proj <- deg2UTM(lonlat[2], lonlat[1]) # Estimate inlet node UTM zone
  # Convert UTM zone to EPSG code. 327xx for south, 326xx for north.
  if(nchar(inlet_node_proj["zone"]) == 1) {
    inlet_node_proj["zone"] <- paste0("0", inlet_node_proj["zone"])
  }
  
  if(inlet_node_proj["h"] == "N") {
    new_epsg <- as.integer(paste0("326",inlet_node_proj["zone"]))
  } else if(inlet_node_proj["h"] == "S") {
    new_epsg <- as.integer(paste0("327",inlet_node_proj["zone"]))
  }
  
  # reference_dot$a_width2 <- reference_dot$apex_width
  # reference_dot$apex_width <- 200L
  
  # Transform inlet node to new coordinates
  reference_dot_UTM <- st_transform(reference_dot, crs = new_epsg)
  # Write inlet nodes, note new filename
  write_sf(reference_dot_UTM, paste0("outlines/",dname,"_reference_dot.shp"),
           delete_layer = T)
  
  ### If you are writing occurrence for first time
  
  
  occurrence_UTM = st_warp(occurrence, crs = new_epsg, cellsize = pixres)
  ## Reproject raster onto a new crs ("warp" the grid), save the UTM for future
  
  write_stars(occurrence_UTM, fn_oc_UTM,
              options = c("COMPRESS=LZW", "TFW=NO"))
  rm(occurrence, new_epsg, lonlat, inlet_node_proj, reference_dot)
} else {
  occurrence_UTM <- read_stars(fn_oc_UTM, proxy = F)
}

# get NA ind of occurrence_UTM, we'll need this later for
# OAM
NA_ind <- which(is.na(as.matrix.stars(occurrence_UTM)))
# Get the occurrence data into R array format 
# stars does something with javascript behind the scenes 
# so we use a wrapper
watermask <- as.matrix.stars(occurrence_UTM)


complex_outer <- T
# set certain values to water, use raster's faster "mask" fcn for this. Future
# scripts should use terra. Most of that is literally going back and forth between 
# matrix, stars, and raster representation of data
if(complex_outer) {
  outer_rm <- read_sf("outlines/outer_rm.shp")
  watermask_UTM <- st_as_stars(watermask)
  st_dimensions(watermask_UTM) <- st_dimensions(occurrence_UTM)
  watermask_rr <- as(watermask_UTM, "Raster")
  watermask_rr <- mask(watermask_rr, outer_rm, updatevalue = 100, inverse = T)
  watermask <- t(as.matrix(watermask_rr))
  watermask_UTM <- st_as_stars(watermask_rr)
  rm(watermask_rr)
  gc()
}

# Threshold used to create watermask (0 to 100% )
tht <- 50

# Threshold the mask and set val > tht to 1
watermask <- bin_thresh(watermask, tht)
mode(watermask) <- 'integer'

if(length(NA_ind) >= 1){
  watermask[NA_ind] <- 0
}

# Only need largest connected component to actually run OAM,
# lakes/interior disconnected channels are computational wastes
channelmask <- keep_largest(watermask)
channelmask <- fill_holes(channelmask, 5) # fill islands less than 5 pixels
# retain binariness of channel mask
channelmask[channelmask>=1] <- 1

# Repeat this after doing keep_largest in case something weird happens
if(complex_outer) {
  outer_rm <- read_sf("outlines/outer_rm.shp")
  channelmask_UTM <- st_as_stars(channelmask)
  st_dimensions(channelmask_UTM) <- st_dimensions(occurrence_UTM)
  channelmask_rr <- as(channelmask_UTM, "Raster")
  channelmask_rr <- mask(channelmask_rr, outer_rm, updatevalue = 1, inverse = T)
  channelmask <- t(as.matrix(channelmask_rr))
  channelmask_UTM <- st_as_stars(channelmask_rr)
  rm(channelmask_UTM,channelmask_rr)
}

# Set values to land - erase water features. This is distinct from "channel_rm" 
# below which is only for preprocessing, this is to actually convert 
manual_cut <- F
if(manual_cut) {
  manual_remove <- st_read("outlines/manual_remove.shp")
  channelmask_UTM <- st_as_stars(channelmask)
  st_dimensions(channelmask_UTM) <- st_dimensions(occurrence_UTM)
  channelmask_UTMr <- as(channelmask_UTM, "Raster")
  channelmask_UTMr <- mask(channelmask_UTMr,manual_remove, updatevalue = 0, inverse = T)
  channelmask <- t(as.matrix(channelmask_UTMr))
  channelmask <- keep_largest(channelmask)
  channelmask[channelmask>0] <- 1
  rm(channelmask_UTMr, channelmask_UTM)
}

# At this stage if you want to just check what the channelmask looks like
# use EBImages javascript based viewer
EBImage::display(channelmask)

# When projecting from EPSG:4326 (WGS84) to a local grid you get "NA" values
# These should be set to 1 for input or zero for saving if want to use it
# for network analysis.
# use = 0 if using custom mask, = 1 if using GSW occurrence
channelmask[NA_ind] <- 0


channelmask_UTM <- st_as_stars(channelmask)
st_dimensions(channelmask_UTM) <- st_dimensions(occurrence_UTM)


watermask_UTM <- st_as_stars(watermask)
st_dimensions(watermask_UTM) <- st_dimensions(occurrence_UTM)

# Write water and channelmasks for reference
fn <- paste0("./thresh",tht,"_channels.tif")
write_stars(channelmask_UTM, fn,
            options = c("COMPRESS=LZW", "TFW=NO"),
            progress = T,
            type = 'Byte')
fn <- paste0("./thresh",tht,"_watermask.tif")
write_stars(watermask_UTM, fn,
            options = c("COMPRESS=LZW", "TFW=NO"),
            progress = T,
            type = 'Byte')

rm(watermask_UTM, channelmask_UTM, watermask)

# For the base GSW occurrence product, they do not compute
# water more than 5 km away from land. Smart! But we need to remove those
# "blind" pixels
bufferissue <- F

##### #####
if(bufferissue) {
  # Identify islands
  channelmask_inv <- (!channelmask)*1
  islands <- label_wrapper(channelmask_inv) 
  
  mode(islands) <- 'integer'
  # Find island which has minimum distance from its shoreline to the cardinal direction
  
  cbuff = card_select(islands, 9)
  
  testmat <- matrix(1, nrow = nrow(islands), ncol = ncol(islands))
  testmat[cbuff[1], cbuff[2]] <- 0
  testmat <- as.matrix(imager::distance_transform(imager::as.cimg(testmat), 0)[, , , 1])
  island_ind <- split_objects(islands)
  to_rm <- sapply(island_ind, function(x) min(testmat[x])) |>
    which.min() |>
    names()
  
  
  islands[unlist(island_ind[to_rm])] <- 0
  islands[islands > 0] <- 1
  channelmask <- (!islands)*1
  channelmask_r[[1]] <- channelmask
  
  display(channelmask)
  rm(channelmask_inv, islands, testmat, cx, cy)
  gc()
  
}

display(channelmask)

# In some cases deltas are way bigger than the relevant shoreline, in those cases
# clip to the shoreline.
toobig_MARK <- F
if(toobig_MARK) {
  toobig <- read_sf("outlines/toobig.shp")
  channelmask_r[[1]][NA_ind] <- 1
  channelmask_rr <- as(channelmask_r, "Raster")
  crs(channelmask_rr) <- st_crs(channelmask_r)$proj4string
  channelmask_rr <- crop(channelmask_rr, toobig)
  NA_ind2 <- which(is.na(t(as.matrix(channelmask_rr))))
  channelmask_r <- st_as_stars(channelmask_rr)
  channelmask <- as.matrix.stars(channelmask_r)
  
  channelmask[NA_ind2] <- NA # new NA ind values to use later on
  display(channelmask)
}


channel_rm <- read_sf("outlines/channel_rm.shp")
# channelmask[NA_ind] <- 1 # 1 or NA depending on base or custom occurrence
channelmask_r[[1]] <- channelmask
channelmask_rr <- as(channelmask_r, "Raster")
channelmask_rr <- mask(channelmask_rr, channel_rm, updatevalue = 0, inverse = T)
channelmask_r <- st_as_stars(channelmask_rr)
channelmask <- as.matrix.stars(channelmask_r)
rm(channelmask_rr)
# channelmask[NA_ind] <- 1 # 1 or NA DOUBLE CHECK. set as NA_ind2 if using toobig
st_crs(channelmask_r) <- st_crs(channel_rm)
display(channelmask)


# For couple of deltas even faster to split the raster into subsections
# Because the OAM is sensitive to what is "seen" nearby the cuts are manually defined

# 1:4273, 4273:8826, 8826:nrow(channelmask_r)
## GBM cuts:
# cuts <- list(a = 1:4273,
#              b = 4273:8826,
#              c = 8826:nrow(channelmask_r))

#Lena
# cuts <- list(a = 1:5229,
#              b = 5229:nrow(channelmask_r)) 
#Yana to really speed it up
# cuts <- list(a = 1:1769,
#              b = 1769:nrow(channelmask_r))
#Niger to speed it up
# cuts <- list(a = 1:3997,
#              b = 3997:nrow(channelmask_r))


# for(k in 1:length(cuts)) {
#   channelmask_test <- channelmask_r[, cuts[[k]], ]
#   channelmask_p <- matrix(as.matrix(channelmask_test[[1]]), nrow = dim(channelmask_test)[1], ncol = dim(channelmask_test)[2])
#   
#   ptm <- proc.time()
#   OAM_map <- ROAM(channelmask_p, precision = 720, save_im = T,
#        fn_r = paste0("OAM_map_", k,".tif"), no_cores = 3,
#        parallel = 2)
#   proc.time() - ptm
#   rm(channelmask_test, channelmask_p, OAM_map)
#   gc()
#   
# }

# Run ROAM-style OAM
ptm <- proc.time()
OAM_map <- ROAM(channelmask, precision = 720, save_im = T,
               fn_r = "OAM_map.tif", no_cores = 3,
               parallel = 2)
# This function will print the number of query points it's evaluating
# in the present iteration (doesn't count # of iterations)
proc.time() - ptm
gc()
beepr::beep(4)


# Extract shoreline in 5 degree critical angle increments, then save as geojson
# Check what cardinal direction to use
theta_list <- seq(10, 90, by = 5)
shoreline_list <- lapply(theta_list, \(th) extract_shoreline(OAM_map, theta = th, card = 5, NA_buff = NA_ind))

shoreline_list2 <- do.call(rbind, shoreline_list)
shoreline_list2$A1 <- theta_list
colnames(shoreline_list2)[colnames(shoreline_list2) == "A1"] <- "theta_c"
shoreline_list2$length <- st_length(shoreline_list2)

# MUST delete underlying json, json's don't overwrite. But 1 file versus 6.
write_sf(shoreline_list2, dsn =  paste0(getwd(), "/outlines/shoreline_variability_walker.json"),
         driver = "GEOJSON", delete_layer = T)
beepr::beep(4)