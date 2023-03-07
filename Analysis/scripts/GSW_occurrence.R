### Download GSW 1984-201X occurrence data

library(reticulate)
library(rgee)
library(sf)
library(geojsonio)


# ee_install(py_env = "rgee")
use_condaenv("rgee", required = T)

ee_Initialize()

gm <- import("geemap")

occur = ee$Image("JRC/GSW1_2/GlobalSurfaceWater")$select(0); # occurrence info
print(occur)

Map$setCenter(-110, 65)

colz <- hcl.colors(100, palette = 'Earth')
visParams = list(
  min = 0.0,
  max = 100,
  bands = "occurrence",
  palette = colz
)

Map$addLayer(
  eeObject = occur,
  visParams = visParams,
  name = "GSW Occurrence"
)

proj <- occur$projection()$getInfo()$crs
scale <- occur$projection()$nominalScale()$getInfo()


## Load up saved outlines (later can add more/less to this crazy shapefile)
doutlines <- read_sf("~/deltas/shoreline_morphometric_space/analysis/globaldata/delta_outlines/delta_outlines.shp")

# Convert to GEE
outlines <- sf_as_ee(doutlines)

Map$setCenter(-110, 65) # Longitude, latitude of where you'll set the GEE map to view

## Use this to plot the image. Visualization is kind of broken
# google hcl.colors to get list of palettes
colz <- hcl.colors(100, palette = 'Earth')
visParams = list(
  min = 0L,
  max = 100L,
  bands = "occurrence",
  palette = colz
)

# add visual parameters
Map$addLayer(
  eeObject = occur,
  visParams = visParams,
  name = "GSW Occurrence") +
  Map$addLayer(outlines)



# Need a list object to iterate over
outlines_list <- outlines$toList(outlines$size())


for(i in seq_len(nrow(doutlines))) {
  # roi <- outlines$filterMetadata('name', 'equals', dname)
  roi = ee$Feature(outlines_list$get(i-1))
  # print(i)
  idn = as.character(outlines_list$get(i-1)$getInfo()$properties$ID)
  print(idn)
  fold_n = dname; # Folder to save occurrence into in google drive
  
  task_img <- ee_image_to_drive(
    occur,
    description = idn,
    folder = fold_n,
    region = roi$geometry(),
    scale = scale,
    crs = proj,
    maxPixels = 1e13,
    fileNamePrefix = idn,
    timePrefix = F
  )
  task_img$start()
  ## If you want to track a certain task, use.
  ## Dont turn this on normally or it'll
  ## wait to finish first retrieival
  ## before going to the next iterator 
  # ee_monitoring(task_img)
  
}