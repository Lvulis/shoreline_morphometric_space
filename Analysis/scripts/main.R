# Extraction and analysis of multiple spatial-series of delta OAM-shorelines
# Last updated: February 9, 2023. 
# V1.0 Author: lawrence vulis lvulis [a] uci.edu


#### Load libraries ####

library(stars) # raster-handling
library(raster) # also raster-handling... but for certain functions
library(sf) # vector handling
library(EBImage) # image processing
library(parallel) # parallelization
library(compiler) # speed
library(biwavelet) # wavelet
library(plotly)    # SMS figure is done in plotly
library(conicfit) # keep in mind can use this to change circle
library(smoothr) # smoothing 
# This script requires smoothr but it isn't loaded to prevent namespace conflicts.


dirname <- "~/deltas/shoreline_morphometric_space/analysis/"

setwd(dirname)

# load utility functions
source("./scripts/utils.R")

#### Load data ####

# sediment flux data
Q_mat = read.csv("sediment_fluxes.csv", header = TRUE)

dname_l <- rownames(Q_mat) <- Q_mat$Name
Q_mat <- subset(Q_mat, select = -Name)
Q_mat[Q_mat<0] <- NA

sediment_data <- Q_mat/rowSums(Q_mat)
sediment_data <- data.frame(sediment_data)

# process shorelines into time-series, 
# Going to store: raw coordinate series, the filtered coordinate series, 
# and the projected series.
rawdat = smoothdat =  circ_srs = vector(mode = "list", length = length(dname_l))
names(smoothdat) = names(circ_srs) = names(rawdat) <- dname_l

# Constants for shoreline data:
pixres = 30 # Images will be in 30-m resolution from Landsat
write_fresh = F # controls whether shapefiles are written for visualization in a GIS
read_mouths = T # If mouths have already been defined, set to T. 
# If mouths haven't been defined, set to F.
# After first time mouths are defind need to manually check them in a GIS (say Q)
# Some deltas have poor resolution of mouths. Fix those manually.
mouth_manual <- c("Nadym", "Mississippi", "Don", "Dnieper",
                  "Apalachicola", "ChaoPhraya", "Atchafayla",
                  "Kobuk", "Hong", "Yukon", "Selenga")

# Use for loop because getting caught in an apply environment is dangerous here
for(i in seq_along(dname_l)) {
  dname = dname_l[i]
  setwd(paste0(dirname, dname))
  
  reference_dot <- st_read(paste0("./outlines/",dname,"_inlet_nodes.shp"))
  
  
  if(file.exists("./outlines/shoreline_variability_walker_rayshader3.json")) {
    shoreline_L2 <- read_sf("./outlines/shoreline_variability_walker_rayshader3.json")
  } else {
    shoreline_L2 <- read_sf("./outlines/shoreline_variability_walker.json")
  }
  
  
  # shoreline_L2 <- read_sf("./outlines/shoreline_variability_walker.json")
  # Some of them ended up unordered for some reason. undo that.
  shoreline_L2 <- shoreline_L2[order(shoreline_L2$theta_c), ]
  
  # Use the 45 degree contour
  theta_select <- 8
  
  # The base coordinate series is sampled on pixel centers which can be pixres or pixres*sqrt(2) apart. 
  # Down sample and resample to 30-m with marginal changes. 
  densified <- smoothr::smooth(shoreline_L2[theta_select, ], method = 'densify', max_distance = pixres / 4)
  xy_orig <- st_coordinates(densified)
  coordinate_series <-  simplify_curve(xy_orig, pixres)[, 1:3] # This is like what smoothr does but custom
  
  # Store the raw data series
  rawdat[[dname]] <- coordinate_series
  
  # Create linestring object from the coordinates
  temp_line <- st_linestring(coordinate_series) |>
    st_sfc() |>
    st_sf() |>
    st_set_crs(st_crs(shoreline_L2)) |>
    st_zm()
  
  
  #Convert to data.frame for next part of adding in mouths
  coordinate_series <- as.data.frame(coordinate_series)
  coordinate_series$RN <- 1:nrow(coordinate_series)
  
  
  
  if(read_mouths) {
    ### Just read them! 
    new_SL = st_read("./outlines/OAM_mouth_inanalysis_90star.shp")
    new_SL <- new_SL[!st_is_empty(new_SL), ]
    new_SL <- st_cast(new_SL, "LINESTRING")
  } else {
    ### Define them! Buffer between the 45\60 and 45\90 sets
    ### Any mouth pieces which intersect both stay, and take the 45\90. 
    # This helps decrease amount of manual clipping but is not a *necessary* addition. 
    SL_90 = st_buffer(tail(shoreline_L2[, 'id'], 1), dist = pixres)
    mouth_90 = st_difference(temp_line, SL_90) |>
      rename_geometry("geometry") |>
      st_cast("LINESTRING")
    
    SL_60 = st_buffer(shoreline_L2[theta_select + 3, 'id'], dist = pixres)
    mouth_60 = st_difference(temp_line, SL_60) |>
      rename_geometry("geometry") |>
      st_cast("LINESTRING")
    
    mouth_90$id <- seq_len(nrow(mouth_90))

    
    test9060 = st_intersects(mouth_90, mouth_60)
    mouth_90 <- mouth_90[lengths(test9060)>0, ] 

    if(dname %in% mouth_manual) {
      mouthends = st_read("./outlines/mouth_fix.shp")
      mouthends_coord = as.matrix(st_coordinates(mouthends))
      new_pieces = do.call(rbind, lapply(unique(mouthends$id), \(idP) {
        idP_iter = which(mouthends$id==idP)
        mouthids = apply(mouthends_coord[idP_iter, ], 1, \(x) {
          sweep(coordinate_series[, 1:2], 2, x[1:2], "-")^2 |>
            rowSums() |>
            sqrt() |>
            which.min()
        })
        
        
        new_pc <- st_linestring(as.matrix(coordinate_series[mouthids[1]:mouthids[2], 1:2])) |>
          st_sfc() |>
          st_sf() |>
          st_set_crs(st_crs(shoreline_L2)) |>
          rename_geometry("geometry")
        new_pc
      } ))
      
      new_pieces$id <- (1:nrow(new_pieces)+max(mouth_90$id))
      
      mouth_90 = rbind(mouth_90, new_pieces)
      covered = lengths(st_covered_by(mouth_90))
      mouth_90 <- mouth_90[covered<2, ]
    }
    
    
    st_write(mouth_90, "./outlines/OAM_mouth_inanalysis_90star.shp",
             delete_layer = T)
    st_write(mouth_60, "./outlines/OAM_mouth_inanalysis_FULL.shp",
             delete_layer = T)
    
    new_SL = mouth_90
  } 
  
  ## get coordinates of mouths
  mouth_segments <- st_coordinates(new_SL)
  # mouth_segments is a bunch of disconnected segments, just keep coord & then tag them with a "1"
  # Doesn't keep track of number of segments etc., just the location of mouth coordinates 
  mouth_segments <- cbind(mouth_segments[, 1:2], 1)
  colnames(mouth_segments)[3] <- 'M' # Call a mouth "1", non-mouth "0"
  
  # Merge the two matrices and reorder by RNto keep coord series order 
  labeled_coordinates_temp = merge(coordinate_series, mouth_segments, all.x=T,sort=F)
  labeled_coordinates_temp$M[is.na(labeled_coordinates_temp$M)] <- 0
  labeled_coordinates_temp <- labeled_coordinates_temp[order(labeled_coordinates_temp$RN), ]
  labeled_coordinates_temp = base::subset(labeled_coordinates_temp, select= -c(RN)) # Drop the order column.
  rm(mouth_segments)
  gc()
  
  # labeled_coordinates_temp now has coordinates and if they're mouths the third column "M" is 1. Need as matrix type for later fcn.
  labeled_coordinates_temp <- as.matrix(labeled_coordinates_temp)
  
  # Smooth coordinate series and resample at sample rate deltas = pixres*2 
  deltas <- pixres*2 
  # Nadaraya-Watson kernel smoothing on x(s) and y(s) independetly
  coordinate_series2 = smoothr::smooth_ksmooth(labeled_coordinates_temp, wrap = F, bandwidth = 6*pixres) |>
    simplify_curve(newstep = deltas)
  colnames(coordinate_series2) = colnames(labeled_coordinates_temp)
  
  # Threshold the mouth column, anywhere its mostly "mouth" from the smoother is mouth, and vice-versa for not mouth
  coordinate_series2[coordinate_series2[, 'M'] < .5, 'M'] <- 0
  coordinate_series2[coordinate_series2[, 'M'] >= .5, 'M'] <- 1

  smoothdat[[i]] <- coordinate_series2
  
  # Write smoothed coordinates to actually do visual inspection & plotting in GIS
  if(write_fresh) {
    coordtest <- coordinate_series2 |>
      st_linestring() |>
      st_sfc() |>
      st_sf() |>
      st_set_crs(st_crs(reference_dot))
    st_write(coordtest, "./outlines/OAM_inanalysis.shp",
             delete_layer = T)
    
  }
  # Get distance along points. Note that a smoothed curve was generated, and 
  # every 60 m along the curve was sampled. Those points may be less
  # than 60 m apart as the crow flies. The distance along the 
  # smoothed curve is assumed to be the correct distance, 
  # rather than computing it from the now downsampled series.
  d_along <- (0:(nrow(coordinate_series2) - 1)) * deltas
  d_along_norm <- d_along / max(d_along)
  
  circfit = circular::lsfit.circle(coordinate_series2[, 1:2])
  centroid_coord_circ <- circfit$coefficients[c(2, 3, 1)]
  
  # Write center of curvature for plotting
  if(write_fresh) {
    geo_centroid_circ <- st_point(circfit$coefficients[2:3]) |>
      st_sfc() |>
      st_sf() |>
      st_set_crs(st_crs(reference_dot))
  st_write(geo_centroid_circ, "./outlines/center_of_curvature.shp",
           delete_layer = T)
  }
  
  circ_srs[[dname]] <- list(cbind(d_along,
                                  d_along_norm,
                                  as.numeric(circfit$radius),
                                  coordinate_series2[, 'M']),
                            circfit[c("coefficients", "angles")])
  rm(d_along, d_along_norm)
  
  setwd("..")
}

### Extract series statistics & information #####

sl_len <- sapply(circ_srs, \(x) {
  x[[1]][nrow(x[[1]]), "d_along"]
})
names(sl_len) <- dname_l
circ_rad <- sapply(circ_srs, \(x) x[[2]]$coefficients[[1]])

#### Plot circle-mapped spatial-series #####

lapply(dname_l, \(dname) {
  
  yl = range(circ_srs[[dname]][[1]][, 3])
  filename = paste0(dirname, "plots/dcircle/", dname, "_45.png")
  png(file = filename, width = 500, height = 300, pointsize = 14)
  par(mar = c(5.1, 5.1, 4.1, 2.1), mfrow = c(1, 1))
  plot(circ_srs[[dname]][[1]][, "d_along"],
       circ_srs[[dname]][[1]][, 3],
       type = 'l',  
       xlim = range(circ_srs[[dname]][[1]][, "d_along"]) + c(7900, 0),
       xlab = bquote("Distance along Shoreline (m)"),
       ylab = bquote(italic(d[c])~"(m)"),
       main = dname,
       cex.lab=1.6,cex.axis=1.3,cex.main=1.6,
       lwd =  2,
       panel.first = grid(lwd = 2, col = rgb(0.05, 0.05, 0.05, 0.5)))#,
  
  IN = circ_srs[[dname]][[1]][, 4] == 0 
  points(circ_srs[[dname]][[1]][IN, 1], circ_srs[[dname]][[1]][IN, 3],
         col = 'black', pch= 16)
  points(circ_srs[[dname]][[1]][!IN, 1], circ_srs[[dname]][[1]][!IN, 3],
         col = 'red', pch= 16)
  rm(IN)
  
  dev.off()
  
})

#### Classify deltas into Concave, Convex, Flat ####

# Convert angular range into ang/2/pi (total circle)
angrange = sapply(circ_srs, \(x) {
  range(x[[2]]$angles)
}) / 2 / pi

# Identify flat deltas
dname_flat = names(which(angrange < 30/360))

# Identify concave deltas (manual). To do automatically could test against casting to center of curvature

dname_concave <- c("Amazon", "Colorado", "Nadym", "Pur",
                   "Dnieper", "Kolyma",  "ChaoPhraya", "Ob",
                   "Yangtze", 'TigrisEuphrates')

# Identify convex deltas (set diff)
dname_convex <- setdiff(dname_l, c(dname_flat, dname_concave))


#### Compute sub-macro scale metrics using wavelets ####

## wavelet transform the series
# The biwavelet package "correction" is an overcorrection and is WRONG
# do NOT use their correction.
# Double correction turns WN into colored noise.
# Thank you to Clement Guilloteau for help with this point.
circ_wt = lapply(seq_along(circ_srs[]), \(i) {
  print(i)
  wavo <- biwavelet::wt(circ_srs[[i]][[1]][, c(1, 3)], 
                        do.sig = F,
                        max.scale = max(circ_srs[[i]][[1]][, c(1)]),
                        mother = 'morlet')
  
  yl <- range(wavo$period)
  yl[2] <- (.35*max(wavo$t)) 
  
  filename = paste0(dirname, "plots/dcircle/Wavelet_", dname_l[i], "_45.png")
  png(file = filename, width = 500, height = 300, pointsize = 14)
  par(mar = c(5.1, 5.1, 4.1, 2.1), mfrow = c(1, 1))
  
  plot(wavo, type = "power.norm",
       main = "Scalogram",
       plot.phase = F,
       fill.cols = hcl.colors(64, palette = 'ag_sunset'),
       xlab = 'Distance Along Shoreline (m)',
       xlim = range(wavo$t),
       ylab = '',
       cex.lab = 1.6, cex.axis = 1.3, cex.main = 1.6,
       las = 1, yaxt = '', ylim = yl)
  
  title(ylab = 'Scale (m)', line = 4, cex.lab = 1.6)
  
  # add mouths as rectangles
  mdiff = diff(circ_srs[[i]][[1]][, 4])
  id_end = which(mdiff<0)
  id_start = which(mdiff>0)
  
  if(id_start[1]>id_end[1]){
    id_start <- c(1, id_start)
  }
  
  sapply(seq_along(id_start), \(j) {
    polygon(x=c(circ_srs[[i]][[1]][id_start[j], 1],
                circ_srs[[i]][[1]][id_end[j], 1],
                circ_srs[[i]][[1]][id_end[j], 1],
                circ_srs[[i]][[1]][id_start[j], 1]),
            y=log2(c(yl[1], yl[1], yl[2], yl[2])),
            col = rgb(.5,.5,.5,.5),
            border = NA)
  })
  
  # Create single rectangle for finescale variance: 
  polygon(x = c(min(wavo$t),
                max(wavo$t),
                max(wavo$t),
                min(wavo$t)),
          y = log2(c(300, 300, 1000, 1000)),
          col = rgb(.5,.5,.5,.5),
          border = NA)
  dev.off()
  
  ### use this once to get a legend manually
  # plot(wavo, type = "power.norm",
  #      plot.phase = F,
  #      fill.cols = hcl.colors(64, palette = 'ag_sunset'),
  #      ylab = '',
  #      cex.lab=1.6,cex.axis=1.3,cex.main=1.6,
  #      las=1, yaxt = 'n', xaxt= 'n', xlab = "", ylim = yl,
  #      legend.horiz=T, plot.cb = T)
  
  return(wavo)
})



#### Wavelet-based PSD ####

## Compute mouth variance ##

mouth_var = sapply(seq_along(circ_wt), \(i) {
  compute_mouth_var(circ_wt[[i]], circ_srs[[i]][[1]][, 4]==0)
})

names(mouth_var) <- dname_l

# Compute finescale variance
circ_wt_int_noCOI <- lapply(circ_wt, wt_cum_psd)

# Define frequency bounds: Change these to make the supplementary figures
lb_fine = 300
ub_fine = 1000
freq_brk <- c(lb_fine, ub_fine) # for smart computation of fine-energy using vector


WT_finenergy_noCOI <- sapply(circ_wt_int_noCOI, \(x) bw_power(x, freq_brk, norm = F))

names(WT_finenergy_noCOI) <- dname_l

# Compute Gini Coefficient over right range of coefficients

WT_finenergy_gini = WT_finenergy_noCOI * sapply(circ_wt_int_noCOI, spectral_gini, freq_brk)

#### Clustering ####


df <- data.frame(FVp_gini = WT_finenergy_gini,
           length = log10(sl_len),
           mouthvar = mouth_var[names(WT_finenergy_gini)])

df$R <- NA
df[dname_convex, "R"] <- 1
df[dname_concave, "R"] <- 3
df[dname_flat, "R"] <- 2

df$RFac = factor(df$R)
levels(df$RFac) <- list("Convex" = "1",
                        "Flat" = "2",
                        "Concave" = "3")
# Rescaled dataframe for clustering (not plotting)

df2 <- df[, c("mouthvar", "FVp_gini")]
df2 <- data.frame(apply(df2, 2, scale))

rownames(df2) <- rownames(df)

df2$R <- factor(df$R)

KM = clustMixType::kproto(df2, k = 5, nstart = 100, iter.max = 1e3)

# Randomly 80% of the deltas and see how well it matches original classificaiton
test_robustness = T
if(test_robustness) {
  LPOCV = lapply(1:100, \(k) {
    
    idkp = sort(sample(1:nrow(df), nrow(df)*.80))
    
    df3 <- df[idkp, c("mouthvar", "FVp_gini")] 
    nm = rownames(df3)  # scale()
    df3 <- data.frame(apply(df3, 2, scale))
    
    rownames(df3) <-nm
    
    df3$R <- factor(df$R[idkp])
    
    KM2 = clustMixType::kproto(df3, k = 5, nstart = 100, iter.max = 1e3)
    plot(df3[, 1:2], pch = df3$R |> as.integer()+15, col = KM2$cluster)
    classd = t(apply(KM2$centers[, 1:2], 1, \(x) {
      sqrt(rowSums((KM$centers[, -which(colnames(KM$centers) == "R")] - rep(as.numeric(x), each = nrow(KM$centers)))^2)) |> 
        which.min()
    }))
    
    newcl = classd[KM2$cluster] 
    
    smp_cl = rep(0, length = nrow(df)) |> setNames(dname_l)
    smp_cl[rownames(df3)] <- newcl
    smp_cl[smp_cl == 0] <- NA
    
    confusion_matrix(as.factor(KM$cluster), as.factor(smp_cl))
    
  })
  
  
  acc = sapply(LPOCV, \(X) {
    if(ncol(X) != nrow(X)) {
      return(NULL )
    } else {  
      sum(diag(X[-nrow(X), -ncol(X)]))/X[nrow(X), ncol(X)]
    }
  })
  acc[!is.null(acc)] |> unlist()
  
}

df[rownames(df2), "cluster"] = KM$cluster


# Relabel numbered clusters using reference deltas
clA = which(df$cluster==df$cluster[which(rownames(df)=='Fly')])
clB = which(df$cluster==df$cluster[which(rownames(df)=='Indus')])
clC = which(df$cluster==df$cluster[which(rownames(df)=='Wax')])
clD = which(df$cluster==df$cluster[which(rownames(df)=='Ebro')])
clE = which(df$cluster==df$cluster[which(rownames(df)=='Waipaoa')])
df$cluster[clE] <- "Wave"
df$cluster[clA] <- "Tide"
df$cluster[clB] <- "River-Tide"
df$cluster[clC] <- "River"
df$cluster[clD] <- "River-Wave"
df$cluster <- factor(df$cluster)




#### Plot shoreline morphometric space (Figure 3)####
# Uses plotly

flab <- list(
  family = "Arial",
  size = 20,
  color = "black"
)

ftick <- list(
  family = "Arial",
  size = 20,
  color = 'black'
)

tf = list(
  family = "sans serif",
  size = 14,
  color = "black")


pal <- c(rgb(218/255, 37/255, 37/255),
         rgb(64/255, 230/255, 9/255),
         rgb(29/255, 29/255, 226/255),
         rgb(128/255, 128/255, 0/255),
         rgb(200/255, 128/255, 128/255))

pal <- setNames(pal, c("Wave", "River", "Tide", "River-Wave", "River-Tide"))


axpar <- list(
  autotick = T,
  tickfont = ftick,
  mirror = TRUE,
  showline = TRUE,
  titlefont = flab,
  zerolinecolor = '#ffff', 
  zerolinewidth = 2, 
  gridcolor = 'ffff'
)

fig <- plot_ly()
fig <- add_trace(fig,
                 data = df,
                 x = ~mouthvar,
                 y = ~FVp_gini,
                 symbol = ~RFac,
                 symbols = c("circle", "square", "triangle-down"),
                 colors = pal,
                 type = 'scatter', 
                 mode = 'markers',
                 showlegend = T,
                 marker = list(color = 'black',
                               size = 12,
                               line = list(color = 'black'))
)

# Deltas where text label has to be moved
out = c("Parana", "Limpopo", "Indigirka", "Irrawaddy",
        "GBM", "Mahakam", "Po", "Orinoco",
        "TigrisEuphrates", "Krishna", "Klamath",
        "Yangtze", "Mekong", "Niger", "Godavari",
        "Ebro", "Ceyhan", "SaoFrancisco",
        "Huanghe", "Danube", "ParanaRO")

fig <- add_trace(fig,
                 data=as.data.frame(df),
                 x = ~mouthvar,
                 y = ~FVp_gini,
                 symbol = ~RFac,
                 color = ~cluster,
                 symbols = c("circle", "square", "triangle-down"),
                 colors = pal,
                 type = 'scatter', mode = 'markers',
                 hoverinfo = rownames(df),
                 showlegend=F,
                 inherit=F,
                 marker = list(size = 12))

fig <- add_text(fig,
                data = df[!(rownames(df) %in% out), ],
                x = ~mouthvar,
                y = ~FVp_gini,
                text = rownames(df[!(rownames(df) %in% out), ]),
                textposition = 'top right',
                textfont = tf,
                showlegend = F,
                inherit = F)

fig <- add_text(fig,
                data = df[(rownames(df) %in% c( "Indigirka",
                                                "Niger",
                                                "Ebro")), ],
                x = ~mouthvar,
                y = ~FVp_gini,
                text = rownames(df[(rownames(df) %in% c('Indigirka',
                                                        "Niger",
                                                        "Ebro")), ]),
                textposition = 'left',
                textfont = tf,
                showlegend = F,inherit = F)

fig <- add_text(fig,
                data = df[(rownames(df) %in% c( "Irrawaddy" ,"Godavari")), ],
                x = ~mouthvar,
                y = ~FVp_gini,
                text = rownames(df[(rownames(df) %in% c("Irrawaddy", "Godavari")), ]),
                textposition = 'right',
                textfont = tf,
                showlegend = F,inherit = F)

fig <- add_text(fig,
                data = df[(rownames(df) %in% c( "Klamath", "SaoFrancisco",
                                                "GBM")), ],
                x = ~mouthvar,
                y = ~FVp_gini,
                text = rownames(df[(rownames(df) %in% c( "Klamath", "SaoFrancisco",
                                                         "GBM")), ]),
                textposition = 'bottom left',
                textfont = tf,
                showlegend = F,inherit = F)

fig <- add_text(fig,
                data = df[(rownames(df) %in% c("Ceyhan", "Huanghe", "Danube", "Krishna",
                                               "TigrisEuphrates")), ],
                x = ~mouthvar,
                y = ~FVp_gini,
                text = rownames(df[(rownames(df) %in% c("Ceyhan", "Huanghe", "Danube", "Krishna",
                                                        "TigrisEuphrates")), ]),
                textposition = 'top left',
                textfont = tf,
                showlegend = F,inherit = F)


fig <- add_text(fig,
                data = df[(rownames(df) %in% c("Yangtze", "Po", "Orinoco",
                                               "Parana", "ParanaRO", "Limpopo")), ],
                x = ~mouthvar,
                y = ~FVp_gini,
                text = rownames(df[(rownames(df) %in% c("Yangtze", "Po", "Orinoco",
                                                        "Parana", "ParanaRO", "Limpopo")), ]),
                textposition = 'bottom right',
                textfont = tf,
                showlegend = F,inherit = F)

fig <- add_text(fig,
                data = df[(rownames(df) %in% c("Mekong")), ],
                x = ~mouthvar,
                y = ~FVp_gini,
                text = rownames(df[(rownames(df) %in% c("Mekong")), ]),
                textposition = 'top',
                textfont = tf,
                showlegend = F,inherit = F)


fig <- add_text(fig,
                data = df[(rownames(df) %in% c("Mahakam")), ],
                x = ~mouthvar,
                y = ~FVp_gini,
                text = rownames(df[(rownames(df) %in% c("Mahakam")), ]),
                textposition = 'bottom',
                textfont = tf,
                showlegend = F,inherit = F)

fig <- layout(fig,  
              xaxis = c(axpar, title = 'Fraction of Variance Contributed by Mouths'),#, type = 'log'),
              yaxis = c(axpar, title = "Fine-Scale Variance (m<sup>2</sup>)"),
              plot_bgcolor= rgb(235/255, 235/255, 235/255, 1),
              legend = list(x = .05, y = .95,
                            bgcolor = rgb(128/255, 128/255, 128/255, 0.25),
                            bordercolor = rgb(235/255, 235/255, 235/255, 0),
                            title=list(text='<b>         Shape</b>',
                                       side = "top",
                                       font = list(family = 'Arial',
                                                   size = 16,
                                                   color = 'black')),
                            font = list(family = 'Arial',
                                        size = 16),
                            itemsizing = 'constant'),
              title = paste0("<b> Shoreline Morphometric Space"),
              showlegend = T) 

fig

origwd <- getwd()
setwd("../paper/figures")
orca(fig, width = 960, height = 493, format = "png", file = "SMS_classification.png",
     scale = 5)
setwd(origwd)
beepr::beep(4)

#### Plot ternary Galloway diagram (Figure 5) ####

ternaxpar <- list(
  autotick = T,
  tickfont = ftick,
  mirror = TRUE,
  showline = TRUE,
  titlefont = flab)
legpar <- list(font = list(
  family = "arial",
  size = 12,
  color = "#000000")
)
regborder <- list(color = 'black',
                  dash = 'dash', 
                  width = 2)


sediment_data[rownames(df), "morphCluster"] <- (df$cluster)

## Combine sediment_data with disturbed in new data frame
## and use that for arrow plotting
fig <- plot_ly(data=data.frame(sediment_data),
               a = ~River,
               b = ~Wave,
               c = ~Tide,
               color = ~morphCluster,
               text = rownames(sediment_data),
               colors = pal,
               type = 'scatterternary', mode = 'markers',
               marker = list(size = 12,
                             symbol = 'circle'))


fig <- layout(fig,
              title = "",
              ternary = list(
                # sum = 1,
                aaxis = c(ternaxpar, title='River'),
                baxis = c(ternaxpar, title='Wave'),
                caxis = c(ternaxpar, title='Tide')
              ),
              legend = list(title = list(text = "<b>Morphotype</b>",
                                         size = 14,
                                         font = ftick),
                            font = list(family = 'Arial',
                                        size = 20,
                                        color = 'black'))
)

fig


origwd <- getwd()
setwd("../paper/figures")
orca(fig, width = 600, height = 600, format = "png", file = "Galloway_SMSclusters.png",
     scale = 5)
setwd(origwd)
beepr::beep(4)


#### Export data to csv ####

df_export = df[, c("length", "RFac", "mouthvar", "FVp_gini", "cluster")]
df_export$length <- 10^df_export$length
df_export$cluster <- as.character(df_export$cluster)

df_export <- cbind(df_export, Q_mat[rownames(df_export), ])
df_export <- df_export[sort(rownames(df)), ]
df_export$length <- df_export$length/1000

write.csv(df_export, "scripts/shorelinemetrics.csv")



#### cPSD for real signals & WN: Supplementary Figure 1 ####


color_scale_vals <- hcl.colors(64, 'viridis')
color_ranks <- round(rank(sl_len))
names(color_ranks) <- names(sl_len)


filename = paste0(dirname, "plots/dcircle/wt_cumspectra_infrequency.png")
png(file = filename, width = 500, height = 500, pointsize = 14)
par(mar = c(5.1, 5.1, 4.1, 2.1), mfrow = c(1, 1))

plot(-50, -50, xlim = c(1 / ub_fine, 1 / lb_fine), ylim = c(1e-3, 1),
     xlab = bquote('Wavenumber'~"(m" ^ {-1} * ")"),
     ylab = 'Cum. Power ',
     cex.lab = 1.6, cex.axis = 1.6, cex.main = 1.6,
     log = '',
     main = 'Finescale Energy')

sapply(seq_along(circ_wt_int_noCOI)[], \(i) {
  
  kp = which.min(abs(circ_wt_int_noCOI[[i]][, 1] - lb_fine))
  kp2 = which.min(abs(circ_wt_int_noCOI[[i]][, 1] - ub_fine))
  bb = circ_wt_int_noCOI[[i]][kp:kp2, 2]
  bb <- bb - bb[1]
  lines(1 / circ_wt_int_noCOI[[i]][kp:kp2, 1],
        bb / (max(bb)),
        col = color_scale_vals[color_ranks[i]],
        lwd = 3, type = 'l')
  
  
  k1 = 1 / circ_wt_int_noCOI[[i]][kp, 1]
  k2 = 1 / circ_wt_int_noCOI[[i]][kp2, 1]
  mx = 1 / (k2 - k1)
  
  lx = 1 / circ_wt_int_noCOI[[i]][kp:kp2, 1]
  ly = (1 / circ_wt_int_noCOI[[i]][kp:kp2, 1]) * mx + -(k1) * mx
  
  lines(lx, ly,
        col = 'black',
        lwd = 2)
  
})

dev.off()

