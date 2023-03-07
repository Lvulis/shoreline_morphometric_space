### GEEE side downloading of monthly water history
### and creating annual/other periodic occurrence/water mask
# rgee::ee_install_set_pyenv(
#   py_path = "C:/Anaconda3/envs/rgee", 
#   py_env = "rgee"
# )
# Minus 1? 
relab <- function(image) {
  # Turns 2 to 1, 1 to 0, and everything else to NULL
  image2 = image$remap( c(2, 1),
                        c(1, 0),
                        NULL,
                        'water')
  # image2 <- image$updateMask(image$select('water')$gt(1))
  # image2 <- image2$expression('war - 1',
  #                             list('war' = image2$select('water')))
  return(image2)
}
# 
# binThreshOcc <- function(image) {
#   
#   # ee$ImageCollection$filter(ee$Filter$gte('month', mo_start)) %>%
#   
#   thres = image.gte(10).rename('thres')
#   return(thres)
#   
# }

library(reticulate)
reticulate::use_condaenv("rgee", required = T)
library(rgee)
library(sf)

gm <- import("geemap")
rgee::ee_Initialize()

MWH = ee$ImageCollection("JRC/GSW1_2/MonthlyHistory") # full stack

proj <- MWH$first()$projection()$getInfo()$crs
scale <- MWH$first()$projection()$nominalScale()$getInfo()


# Identify delta which you wanna extract
# dname = "Pearl"; ## done on one-by-one basis for this monthly water masks
fold_n = paste0(dname, '_GSW_4326_v1_1'); # Folder to save occurrence into in google drive


## Load saved outlines
doutlines <- read_sf("~/deltas/shoreline_morphometric_space/analysis/globaldata/delta_outlines/delta_outlines.shp")
doutlines$name[which(doutlines$name=='Godai')] <- 'Godavari'
# Convert to GEE
outlines <- sf_as_ee(doutlines)
# #Select and create new featureCollectio
roi <- outlines$filterMetadata('name', 'equals', dname)

### Download monthly water masks from the range specified above
# mo_start = 1
# mo_end = 12
# yr_start = 2015
# yr_end = 2018
# 
# 
# MWH_cut = MWH %>%
#   ee$ImageCollection$filter(ee$Filter$gte('month', mo_start)) %>%
#   ee$ImageCollection$filter(ee$Filter$lte('month', mo_end)) %>%
#   ee$ImageCollection$filter(ee$Filter$gte('year', yr_start)) %>%
#   ee$ImageCollection$filter(ee$Filter$lte('year', yr_end))
# 
# 
# 
# MWH_list <- MWH_cut$toList(MWH_cut$size())
# Nimg <- MWH_cut$size()$getInfo()
# for(i in seq_len(Nimg)) {
#   print(i)
#   idn = ee$Image(MWH_list$get(i-1))$getInfo()$id
#   idn <- substr(idn, 27, 33)
#   task_img <- ee_image_to_drive(
#     ee$Image(MWH_list$get(i-1)),
#     description = idn,
#     folder = fold_n,
#     region = roi$geometry(),
#     scale = scale,
#     crs = proj,
#     maxPixels = 1e13,
#     fileNamePrefix = idn,
#     timePrefix = F
#   )
#   task_img$start()
# 
# }



### Generating annual occurrence map or water mask server-side

# yr <- 2017
for(yr in 2010:2020) {
# for(mo in 1:12) {
  print(yr)
  MWH_cut = MWH %>%
    ee$ImageCollection$filter(ee$Filter$eq('year', yr))
  
  
  # Apply relab to the entire image set
  MWH_cut2 <- MWH_cut$map(relab)
  # Apply pixel-based average
  new_occ <- MWH_cut2$reduce(ee$Reducer$mean())$reproject(proj, scale = scale)
  
  # new_occ <- MWH_cut$reduce(ee$Reducer$mean()) 
  new_occ <- new_occ$expression('occ * 100',
                                list('occ' = new_occ$select('remapped_mean'))) %>%
    ee$Image$round()
  
  # thres = new_occ$gte(50)$rename('thres')
  idn <- paste0("occurrence_", yr)
  # idn <- paste0("occurrence_", month.name[mo])
  # idn <- paste0("watermask_", yr)
  
  task_img <- ee_image_to_drive(
    new_occ,
    # thres,
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
  # ee_monitoring(task_img)
}