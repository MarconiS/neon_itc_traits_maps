clean_spectra <- function(brick){
  # filter for no data 
  mask = brick > 10000
  brick[mask] <-  NA
  mask = brick == -9999
  brick[mask] <-  NA
  
  #filter for greennes and shadows
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
  
  #filter for known artifacts
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
  
  # # save pixel positions
  good_pix = !is.na(brick)
  good_pix = (apply(good_pix, 1, all))
  # 
  brick = brick[complete.cases(brick),]
  return(list(refl = brick, good_pix = good_pix))
}

