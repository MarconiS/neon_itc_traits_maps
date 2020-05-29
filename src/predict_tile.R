predict_tile <- function(tile_path, mods, weights){
  tile_path = "../../Chapter4/NEON_D01_HARV_DP3_726000_4699000_reflectance.tif"
  tile_meta = strsplit(tile_path, split = "_")[[1]]
  
  #name of traits predicted 
  tr_nm = c("N", "C", "Crt", "CAC", "CBC", "dryM", "LMA", "lign", "cell")
  brick = raster::brick(tile_path) 
  tile_dim = dim(brick)
  brick = raster::as.matrix(brick)
  colnames(brick) = paste("band", 1:369, sep="_")
  
  #clean tile  
  mask = brick > 10000
  brick[mask] <-  NA
  mask = brick == -9999
  brick[mask] <-  NA
  
  ndvi <- (brick[,"band_90"]- brick[,"band_58"])/(brick[,"band_58"] + brick[,"band_90"]) <0.3
  nir860 <- (brick[,"band_96"] + brick[,"band_97"])/20000 < 0.1
  mask = as.logical(ndvi | nir860)
  brick[mask,] = NA
  rm(mask, ndvi, nir860)
  brick = brick[,15:365]
  normMat <- sqrt(apply(brick^2,FUN=sum,MAR=1, na.rm=TRUE))
  normMat <- matrix(data=rep(normMat,ncol(brick)),ncol=ncol(brick))
  brick=brick/normMat
  rm(normMat)
  
  cnd = (brick[,"band_312"] > 0.025)  
  #idx <- which(apply(cnd, 1, any))
  brick[cnd,] = NA
  
  cnd = (brick[,24:45] > 0.04)
  idx <- (apply(cnd, 1, any))
  if(length(idx) !=0){
    brick[idx,] = NA
  }
  cnd = (brick[,100:200] > 0.15)
  idx <- (apply(cnd, 1, any))
  if(length(idx) !=0){
    brick[idx,] = NA
  }
  rm(cnd,idx)
  
  # save pixel positions
  good_pix = !is.na(brick)
  good_pix = (apply(good_pix, 1, all))
  
  brick = brick[complete.cases(brick),]
  
  # reduce into 1st PC of the KLD 
  cfc_pt = "./outdir/hdr_cfc.rds"
  cfc_reduced_dims = readRDS(cfc_pt)
  brick = cbind.data.frame(cfc_reduced_dims$bands_grouping, t(brick))
  colnames(brick)[1] = "kld_array"
  kld_refl = list()
  #loop through groups and create a PCA for each
  for(gp in unique(cfc_reduced_dims$bands_grouping)){
    pcx = brick %>% dplyr::filter(kld_array == gp) #%>% t %>% prcomp()
    pcxgrp = predict(cfc_reduced_dims$pcas[[gp]], newdata = t(pcx))
    kld_refl[[gp]] = pcxgrp[,1]
  }
  brick = do.call(cbind.data.frame, kld_refl)
  rm(kld_refl, pcx, pcxgrp)
  colnames(brick) = paste("kd", 1:ncol(brick), sep="_")
  
  # add accessory features
  brick = brick[-1,]
  brick["siteID"] = tile_meta[3]
  
  chm_path = list.files('./indir/CHM/', pattern = paste(tile_meta[5], tile_meta[6], sep="_"), full.names = T)
  elev_path = list.files('./indir/DTM/', pattern = paste(tile_meta[5], tile_meta[6], sep="_"), full.names = T)
  
  elevatn = raster::brick(elev_path) 
  elevatn = raster::as.data.frame(elevatn)
  elevatn = chm[good_pix,]
  brick = cbind.data.frame(elevatn, brick)
  
  #scale features
  scaling_fct = readRDS("./indir/scaling_fct.rds")
  for(x in 1:nrow(scaling_fct)){
    tmp = scale(brick[[x]], center = scaling_fct[x,"center"], scale = scaling_fct[x,"scale"])
    brick[[x]] = tmp
  }
    
  #predict tile
  ptile = predict(mods$mod, newdata = brick,
                  nsamples =200, robust = T)
  saveRDS(ptile, paste(outdir, paste(c(tile_meta[c(3:6)], "maps.rds"), collapse ="_"), sep="/"))
  # make it a raster back (include only mu and sd)
  brick = raster::brick(chm_path) 
  lyr = (matrix(NA, 1000000,9))
  lyr[good_pix,] = ptile[,1,]
  dim(lyr) = c(1000, 1000,9)
  lyr = raster::brick(lyr, xmn=brick@extent[1], xmx=brick@extent[2], #nl = 9,
       ymn=brick@extent[3], ymx=brick@extent[4], crs=brick@crs, transpose=FALSE)
  names(lyr) = paste(tr_nm, "mu", sep="_")
  brick <- raster::stack(brick, lyr)
  
  #add uncertainty estimates
  lyr = (matrix(NA, 1000000,9))
  lyr[good_pix,] = ptile[,2,]
  dim(lyr) = c(1000, 1000,9)
  lyr = raster::brick(lyr, xmn=brick@extent[1], xmx=brick@extent[2], #nl = 9,
                      ymn=brick@extent[3], ymx=brick@extent[4], crs=brick@crs, transpose=FALSE)
  names(lyr) = paste(tr_nm, "sd", sep="_")
  brick <- raster::stack(brick, lyr)
  raster::writeRaster(brick, 
                      paste(outdir, "/tiff/", paste(c(tile_meta[c(3:6)], "maps.tiff"), collapse ="_"), sep="/"))
}
