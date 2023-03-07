## Utility functions for shoreline metrics

#### geospatial ####

rename_geometry <- function(g, name){
  # Rename the geometry column of sf `g` into `name``
  # Arguments:
  #  g: sf object
  #  name: new geometry column name (string)
  # Returns:
  #  x: sf object
  current = attr(g, "sf_column")
  names(g)[names(g)==current] = name
  sf::st_geometry(g)=name
  g
}

deg2UTM <- function(lat, lon) {
  # Converts latitude/longitude to UTM coordinates based on Karney's 2011 method.
  # Taken from Veness (c) Chris Veness 2014-2017
  # Args:
  #  lat: Latitude (from 85 S to 85 N )
  #  lon: Longitude ( from -180 to 180)
  # Returns:
  #  out: Named character with UTM zone, N/S, x, y, convergence, and scale
  #
  falseEasting = 500e3
  falseNorthing = 10000e3
  zone = floor((lon+180)/6)+1 # longitudinal zone
  lambda_0 = deg2Rad(((zone-1)*6 - 180 + 3))
  
  phi = deg2Rad(lat) # Latitude from equator in radians
  lambda = deg2Rad(lon) - lambda_0 # longitude from central meridian
  
  # Constants given by the WGS 84 datum
  a = 6378137
  b = 6356752.314245
  f = 1/298.257223563
  
  k0 = 0.9996 # UTM scale on central meridian
  ## easting, northing: Karney 2011 Eq 7-14, 29, 35:
  e = sqrt(f*(2-f)) # eccentricity
  
  n <- sapply(1:6, function(i) (f/(2-f))^i)
  # For iterating
  nseq <- 1:6
  cl = cos(lambda); sl = sin(lambda); tl = tan(lambda) # trig of the lon
  
  tau = tan(phi) # Tan of the latitude
  
  sig = sinh(e*atanh(e*tau/sqrt(1+tau^2)))
  taup = tau*sqrt(1+sig^2) - sig*sqrt(1+tau^2)
  
  zetap = atan2(taup, cl)
  etap = asinh(sl/sqrt(taup^2 + cl^2))
  
  A = a/(1+n[1]) * (1 + 1/4*n[2] + 1/64*n[4] + 1/256*n[6])
  
  alpha =   c( 1/2*n[1] - 2/3*n[2] + 5/16*n[3] + 41/180*n[4] - 127/288*n[5] + 7891/37800*n[6],
               13/48*n[2] -  3/5*n[3] + 557/1440*n[4] +     281/630*n[5] - 1983433/1935360*n[6],
               61/240*n[3] -  103/140*n[4] + 15061/26880*n[5] +   167603/181440*n[6],
               49561/161280*n[4] - 179/168*n[5] + 6601661/7257600*n[6],
               34729/80640*n[5] - 3418889/1995840*n[6],
               212378941/319334400*n[6])
  zeta = zetap + sum(alpha[nseq] * sin(2*nseq*zetap) * cosh(2*nseq*etap))
  
  eta = etap + sum(alpha[nseq] * cos(2*nseq*zetap) * sinh(2*nseq*etap))
  
  x = k0 * A * eta
  y = k0 * A * zeta
  
  # convergence
  pp = 1 + sum(2*nseq*alpha[nseq] * cos(2*(nseq)*zetap) *cosh(2*nseq*etap))
  qp = sum(2*nseq*alpha[nseq] * sin(2*(nseq)*zetap) *sinh(2*nseq*etap))
  
  gammap = atan(taup/sqrt(1+taup^2)*tl)
  gammapp = atan2(qp, pp)
  gam = gammap+gammapp
  
  # scale
  sphi = sin(phi)
  kp = sqrt(1 - e^2*sphi^2) * sqrt(1 + tau^2) / sqrt(taup^2 + cl^2)
  kpp = A / a * sqrt(pp^2 + qp^2)
  
  k = k0 *kp * kpp
  
  # Shift
  x = x + falseEasting # shift relative to false easting
  if (y < 0 ) y = y + falseNorthing # make y in southern hemisphere relative to false northing
  # rounding
  x = as.numeric(format(x, nsmall = 6))
  y = as.numeric(format(y, nsmall = 6))
  convergence = as.numeric(format(rad2Deg(gam), nsmall = 6))
  scale = as.numeric(format(k, nsmall = 12))
  if(lat >= 0) {
    h = "N"
  } else {
    h = "S"
  }
  out <- c(zone, h, x, y, convergence, scale)
  names(out) <- c("zone", "h", "x", "y", "convergence", "scale")
  out
}

deg2Rad <- function(deg) {
  deg*pi/180
}

rad2Deg <- function(rad) {
  rad*180/pi
}

as.matrix.stars <- function(x) {
  # Convert the first layer of a stars object into a matrix
  # Necessary for some kind of javascript error that arises sometimes
  # Arguments:
  #  x: stars object on a regular grid
  #  i: band
  # Returns:
  #  matrix containing the gridded values
  return(matrix(as.matrix(x[[1]]), nrow = dim(x)[1], ncol = dim(x)[2]))
}

card_select <- function(newmap, card) {
  # Obtain the coordinate of in the "card" position along edge of a matrix
  # Arguments:
  #  newmap: n x m matrix
  #  card: cardinal location, see guideline below
  # Returns:
  #  2-vector with cx, cy.
  if(card == 1) {
    cx = 1; cy = 1
  } else if(card == 2) {
    cx = round(nrow(newmap)/4); cy = 1
  } else if(card == 3) {
    cx = round(nrow(newmap)/2); cy = 1
  } else if(card == 4) {
    cx = round(nrow(newmap)*3/4); cy = 1
  } else if(card == 5) {
    cx = nrow(newmap); cy = 1
  } else if(card == 6) {
    cx = nrow(newmap); cy = round(ncol(newmap)/4)
  } else if(card == 7) {
    cx = nrow(newmap); cy = round(ncol(newmap)/2)
  } else if(card == 8) {
    cx = nrow(newmap); cy = round(ncol(newmap)*3/4)
  } else if(card == 9) {
    cx = nrow(newmap); cy = ncol(newmap)
  } else if(card == 10) {
    cx = round(nrow(newmap)*3/4); cy = ncol(newmap)
  } else if(card == 11) {
    cx = round(nrow(newmap)/2); cy = ncol(newmap)
  } else if(card == 12) {
    cx = round(nrow(newmap)/4); cy = ncol(newmap)
  } else if(card == 13) {
    cx = 1; cy = ncol(newmap)
  } else if(card == 14) {
    cx = 1; cy = round(ncol(newmap)/4)
  } else if(card == 15) {
    cx = 1; cy = round(ncol(newmap)/2)
  } else if(card == 16) {
    cx = 1; cy = round(ncol(newmap)*3/4)
  }
  matrix(c(cx, cy), ncol = 2)
} 

simplify_curve <- \(xy, newstep) {
  # For a 2d planar curve, resample to have d(arc length) constant
  # Will linearly interpolate between existing sample points.
  # Designed with overly-dense points in mind, otherwise the linear
  # interpolation will be crappy.
  # If a function along the curve exists, also interpolate its values. 
  # The matrix *MUST* pass the coordinates in R^2 in the first 2 columns
  # Arguments:
  #  xy: nxP matrix of coordinates in 2d space
  #  newstep: new step size (in arc length/ distance along curve)
  # Returns: 
  #  nx2 matrix of new coordinates
  ds <- c(0, cumsum(point_distance(xy[, 1:2, drop = F])))
  # ds2 <- c(0, cumsum(sqrt(rowSums(diff(xy[, 1:2, drop = F])^2))))
  new_ds <- seq(ds[1], ds[length(ds)], by = newstep)
  
  # Linearlly interpolate new points as necessary
  newxy <- t(sapply(new_ds[-1], \(x) {
    # Find 2 nearest points along curve
    smp <- ds-x
    newid <- which.min(abs(smp))
    # Make sure the "shortest distance" as the one above zero to simplify later calcs
    if(smp[newid] < 0) newid <- newid+1L
    # We are always then taking newid as i and doing bwd diff)
    # For a curve r(s) = <x(s), y<s>), will interpolate b/w r(s_o) and r(s_1):
    # with <x(a), x(a)> = <m_x * (s_a - s_o) + x_o, m_y * (s_a - s_o) + y_o>
    # i.e. r(s_a) = m * (s_a - s_o) + r(s_o)
    # where m = <x_1 - x_o, y_1 - y_o>/(s_1-s_o)
    interend <- xy[(newid-1L):newid, , drop = F] # (r(s_0), r(s_1))
    slopes <- matrix(diff(interend)/diff(ds[(newid-1L):newid]), nrow = 1)
    # Note pass -smp because it's giving s_0 - s_a
    slopes * -smp[newid-1L] + interend[1, , drop = F]
    
  }))
  newxy <- rbind(xy[1, , drop = F], newxy)
  newxy 
}
#### image processing ####

bin_thresh <- function(x, thresh) {
  # Binary threshold a numeric/integer vector or matrix
  # Arguments:
  #  x: Matrix
  #  thresh: threshold value
  # Returns:
  #  x: Binary matrix
  x[x < thresh] <- 0L
  x[x >= thresh] <- 1L
  mode(x) <- 'integer'
  x
}

split_objects = function(x) {
  # From EBImage: get R indices of each object in x matrix
  # Arguments:
  #  x: Integer valued labelled image
  # Returns:
  #  list w/ each entry corresponding to indices of a unique label in x
  z = which(as.integer(x)>0L)
  split(z, x[z])
}

keep_largest <- function(img) {
  # Keep only largest element of an image. Relies on python
  # Arguments:
  #  img: integer-valued binary matrix
  # Returns:
  #  Same sized matrix only containing largest cluster (label id still there)
  cncmp <- label_wrapper(img, conn = T)
  
  mode(cncmp) <- 'integer'
  xid <- split_objects(cncmp)
  ox <- lengths(xid)
  
  if(length(ox) > 1) {
    tokp <- names(which.max(ox))
    cncmp[unlist(xid[setdiff(names(ox), tokp)])] <- 0
  }
  return(bin_thresh(cncmp, 1))
}

label_wrapper <- function(img, conn = T) {
  # Emulate EBImage's bwlabel behavior, s.t. pixels == 0 should appear as 0 in the labelled scene!
  # Arguments:
  #  img: image to be labelled
  #  conn: whether to use low or high connectivity in imager::label
  # Returns:
  #  labelled image
  
  cimgd <- imager::as.cimg(img) > .5
  lab <- imager::label(cimgd, conn) + 1 # |>
  lab <- lab * imager::as.cimg(cimgd)
  
  return(as.matrix(lab[, , , 1]))
}

fill_holes <- function(scene, maxholesize = 1) {
  # Translated from Jon Schwenk, RivGraph.
  # Fills holes of size in binary mask of 0 < size <= maxholesize
  
  cncmp <- label_wrapper((!scene)*1L, conn = F)
  ox <- lengths(split_objects(cncmp))
  
  cncmp[unlist(ox[names(which(ox <= maxholesize))])] <- 0
  cncmp[cncmp > 0] <- 1L
  scene <- (!cncmp)*1L
  scene
}

#### 2-D curve to 1-D series transformations ####
## Construct curvature using equation of JS et al. 2015
## these are setup so for a series from 1 to N, 
## 1 and N will be discarded

## x_diff are helper functions for computeCurvature
bw_diff <- function(M) {
  M[2:length(M)] - M[1:(length(M)-1)]
}
cen_diff <- function(M) {
  M[3:length(M)] - M[1:(length(M)-2)]
}
fw_diff <- function(M) {
  M[3:length(M)] - M[2:(length(M)-1)]
}

compute_curvature <- function(xy, method = 3) {
  # Compute curvature of coordinate pairs {xy} using methods detailed in Schwenk et al., 2015
  # Arguments:
  #  xy: numeric, n x 2 set of coordinates
  #  method: Method to use (see Schwenk et al., 2015)
  # Returns:
  #  Curvature spatial/time-series
  
  if (method == 1) {
    dxy <- diff(xy)
    ds <- sqrt(rowSums(dxy^2))
    ddxy <- dxy/ds
    ds2 <- (ds + c(ds[-1], ds[1]))/2    # Change in ddx per unit length
    Cxy <-  diff(ddxy)/ds2[-1]   #
    k <- (ddxy[-1, 2] * Cxy[, 1] - ddxy[-1, 1] * Cxy[, 2])/
      ((ddxy[-1, 1]^2 + ddxy[-1, 2]^2)^(3/2))
  } else if(method == 3) {
    # Compute curvature with equation given in Schwenk et al. 2015
    a <- apply(xy[-nrow(xy), ], 2, bw_diff)
    b <- apply(xy, 2, cen_diff)
    c <- apply(xy, 2, fw_diff)
    
    k = (2 * (a[, 2]*b[, 1] - a[, 1]*b[, 2])) / sqrt((a[, 1]^2 + a[, 2]^2)*(b[, 1]^2 + b[, 2]^2)*(c[, 1]^2 + c[, 2]^2))
    return(k)}
}

# Arc length essentially from pt to pt
# Copied directly from smoothr
point_distance <- function(x) {
  d <- diff(x)^2
  sqrt(rowSums(d))
}

skewness <- function(x, FP = T, na.rm = F) {
  # Compute skewness coefficient with optional fisher-pearson correction to be unbiased
  # Arguments:
  #  x: numeric vector of data
  #  FP: logical whether to use Fisher Pearson corrxn 
  #  na.rm: NA handling
  # Returns: 
  #  scalar skewness value
  if(na.rm) {
    x <- x[!is.na(x)]
  }
  n<-length(x)
  m3 <- (x - mean(x)) ^ 3
  sk <- (sum(m3) / n) / (sd(x) ^ 3)
  if(FP == T) {
    sk <- (n ^ 2 / ((n - 1) * (n - 2))) * sk
  }
  return(sk)
}

sum_stat <- function(x, na.rm = F) {
  # Compute mean, median, and sd of a vector
  # Arguments:
  #  x: numeric vector of data
  #  na.rm: how to handle NA
  # Returns:
  #  4-vector of statistics
  c(mean = mean(x, na.rm = na.rm), 
    median = median(x, na.rm = na.rm), 
    sd = sd(x, na.rm = na.rm),
    skew = skewness(x, na.rm = na.rm))
}

confusion_matrix <- function(classX, classY) {
  # compute confusion matrix for 2 classifications
  levX = levels(classX)
  levY = levels(classY)
  
  out = matrix(NA, nrow = length(levX)+1, ncol = length(levY)+1)
  for(i in seq_along(levX)) {
    for(k in seq_along(levY)) {
      out[i, k] = length(which((classX==levX[i]) & (classY==levY[k])))
    }
  }
  out[, ncol(out)] = rowSums(out,na.rm=T)
  out[nrow(out), ] = colSums(out,na.rm=T)
  
  rownames(out) = c(levX, "Total")
  colnames(out) = c(levY, "Total")
  
  return(out)
}

#### Shoreline metrics//helpers ####
wt_psd = function(wavo, coi = TRUE, mouths = NULL) {
  # Compute power spectral density from wavelet transform
  # Arguments:
  #  wt: biwavelet object
  #  coi: logical, exclude the COI from the computation, default is TRUE.
  #  mouths: numeric, specifying which points in time (space) correspond to mouths
  # Returns:
  #  n x 2 numeric containing power spectral 
  
  # Filter out COI if necessary
  if(coi) {
    coi_mask = sapply(seq_along(wavo$coi), \(kk) {
      wavo$period>wavo$coi[kk]
    })
    wavo$power[coi_mask] <- NA
  }
  if(!is.null(mouths)) {
    wavo$power[, mouths] <- NA
  }
  # Return 
  PW_full = rowMeans(wavo$power, na.rm = T)*wavo$dt
  
  return(cbind(freq = 1/wavo$period, power = PW_full))
  
}

wt_cum_psd = function(wavo, coi = TRUE, mouths = NULL) {
  # Compute cumulative power spectral density from wavelet transform
  # Arguments:
  #  wt: biwavelet object
  #  coi: logical, exclude the COI from the computation, default is TRUE.
  #  mouths: numeric, specifying which points in time (space) correspond don't correspond mouths
  # Returns:
  #  n x 2 numeric containing power spectral 
  psd = wt_psd(wavo, coi = coi, mouths = mouths)
  intPW2 = -pracma::cumtrapz(psd[, 1], psd[, 2])
  return(cbind("period" = 1/psd[, 1], "cPSD" = intPW2))
}

bw_power <- function(psd, brks, norm = T) {
  # Compute power within given bandwidths. Note to make this compatible
  # with outputs from biwavelet, can pass in frequencies just recall
  # that integration must be done on frequencies.
  # Arguments:
  #  psd: n x 2 numeric, contains frequency OR period and associated power. 
  #  brks: numeric, contains frequencies OR period between which to compute power. At least 2 values must be specified
  #  norm: logical, normalize to have unit power. Default is True
  # Returns:
  #  Power within certain range
  if(length(brks) < 2) {
    stop("Need at least 2 frequencies to compute power between")
  }
  
  brks_ind = sapply(brks, \(xx) {
    bp = abs(psd[, 1] - xx)
    which.min(bp)
  })
  
  diff(psd[brks_ind, 2])
}

compute_mouth_var = function(wavo, mouths) {
  # Compute fraction of variance contributed by mouths
  # Arguments:
  #  wavo: biwavelet object
  #  mouths: logical vector, location of NONMOUTHs
  # Returns:
  #  fraction of variance contributed by mouths
  mouthonly = wt_cum_psd(wavo, coi = FALSE, mouths = mouths)
  all = wt_cum_psd(wavo, coi = FALSE)
  
  # wt_cum_psd normalizes power at each frequency by number of obs
  # which is "wrong" if only integrating over mouths! 
  # Both should normalize by N_obs in time, which cancel out,
  # so simply undo the averages
  return(tail(mouthonly[, 2], 1) / tail(all[, 2], 1) * (sum(!mouths) / length(mouths)))
}


spectral_gini <- function(cPSD, brks) {
  # Compute spectral Gini Coefficient between 2 frequencies. 
  # Arguments:
  #  cPSD: nx2 numeric containing cumulative power spectral density. Defined over PERIODS in first column for convenience
  #  brks: 2-vector containing shorter and longer wavenumber respectively
  # Returns:
  #  Gini coefficient between brks. 
  # Note this function is hardcoded to be for only a single bandwidth
  lb_fine = brks[1]
  ub_fine = brks[2]
  
  kp = which.min(abs(cPSD[, 1]-lb_fine))
  kp2 = which.min(abs(cPSD[, 1]-ub_fine))
  
  # Pull out the cPSD between the reference scales into separate objects
  xb = cPSD[kp:kp2, 1] #scale
  bb = cPSD[kp:kp2, 2] #cPSD
  bb <- bb - bb[1] # reset to zero
  
  # integrate entire bandwidth with respect to frequency
  totint = -pracma::trapz(1 / xb, bb)
  # divide integral by area under WN (line) for same frequency band
  # use 1 - to get complimentary area (gini)
  return(1 - (totint / -(((1 / ub_fine) - (1 / lb_fine)) * diff(cPSD[c(kp,kp2), 2]) / 2)))
  
}
