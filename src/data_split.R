library(tidyverse)

plot_spectra<-function(plt_dat){
  #plot reflectances
  plot_data <- plt_dat %>% 
    dplyr::select(contains("band_")) 
  plot_data <- plot_data %>% dplyr::select(-one_of("individualID")) %>%
    t %>%
    data.frame
  colnames(plot_data) = unlist(plt_dat[1]) # the first row will be the header
  plot_data <- data.frame(bnd = 1:dim(plot_data)[1], plot_data)
  ggdat <- tidyr::gather(plot_data, treeID, Reflectance,-bnd)
  
  return(ggplot(ggdat, aes(x = bnd, y = Reflectance)) + 
           geom_line(aes(color = factor(treeID), alpha= 1), size = 0.2) +
           theme_bw()+
           theme(legend.position="none"))
  
}


optical.filter<- function(dat){
  ndvi <- (dat$band_90- dat$band_58)/(dat$band_58 + dat$band_90) <0.3
  nir860 <- (dat$band_96 + dat$band_97)/20000 < 0.1
  naval = as.logical(ndvi | nir860)
  dat[naval,] = NA
  return(dat)
}

mean.no.na <- function(x){mean(x, na.rm = T)}

hiper_features <- function(dat, normalization = "norm2", type = "spectra"){
  if(type=="spectra"){
    dat <- dat[colnames(dat) %in%
                 paste("band_", seq(1,369), sep="")] %>%
      as.matrix
  }else{
    dat <- as.matrix(dat)
  }
  if(normalization == "norm2"){
    normMat <- sqrt(apply(dat^2,FUN=sum,MAR=1, na.rm=TRUE))
    normMat <- matrix(data=rep(normMat,ncol(dat)),ncol=ncol(dat))
    normMat=dat/normMat
  }else if(normalization == "max_min"){
    min_x = apply(dat, 1, min)
    max_x <- apply(dat, 1, max)
    normMat = dat
    for(jj in 1: dim(normMat)[1]){
      #normMat[jj,] <- (normMat[jj,]-min_x[jj])/(normMat[jj,]-max_x[jj])
      normMat[jj,] <- (normMat[jj,]-min_x[jj])/(max_x[jj]-min_x[jj])
      
    }
    # normMat <- apply(dat,FUN=function(x, na.rm){(x-min_x)/(max_x-min_x)},MAR=1, na.rm=TRUE)
    #normMat <- t(normMat)
  }else if(normalization == "max_min_bnd"){
    scaled <- as.data.frame(scale(dat, center = mins, scale = maxs - mins))
  }else{
    normMat = dat
  }
  return(normMat)
}
#"./ccbid/support_files/Datasets/itc_reflectance.csv"
#/Users/sergiomarconi/Documents/GitHub/hiPyRneon/outdir/spectra/itc2px_reflectance.csv
get_full_dataset <- function(spectra){
  #spectra = readr::read_csv("/Users/sergiomarconi/Documents/GitHub/hiPyRneon/outdir/spectra/itc2px_reflectance.csv")
  spectra = spectra %>% dplyr::select(matches("individualID|band"))
  foo=spectra[-1]
  mask = foo > 10000
  foo[mask] <-  NA
  spectra[-1]=foo
  spectra = dplyr::na_if(spectra, -9999)
  spectra = spectra[complete.cases(spectra),]
  spectra=optical.filter(spectra)
  spectra = spectra[complete.cases(spectra),]
  spectra = spectra[c(1,15:365)]
  spectra[-1] = hiper_features(spectra[-1])
  #spectra[-1] = scale(spectra[-1])
  #vst = sf::read_sf("./indir/traits_polygon_dataset.shp") %>%
  # dplyr::select(indvdID, domanID, siteID, taxonID, scntfcN, elevatn, nlcdCls, cllctDt, ntrgnPr, 
  #               crbnPrc, CNratio , extrcCC, extrCAC, extrCBC, dryMass, dryMssF, 
  #               lfMssPA, leafAre, lgnnPrc, clllsPr, plntStt)
  
  vst = sf::read_sf("./indir/NEON_polygon_dataset.shp") %>%
    dplyr::select(indvdID, domanID, siteID, taxonID, scntfcN)  %>% filter(siteID == "HARV")
  colnames(vst)[1] = "individualID"
  #vst$cllctDt = lubridate::month(vst$cllctDt)
  vegetation_dataset = inner_join(vst, spectra) %>% unique %>% data.frame
  #vegetation_dataset = vegetation_dataset[complete.cases(vegetation_dataset),]
  #readr::write_csv(vegetation_dataset, "~/Documents/Data/vegetation_dataset.csv")
  return(vegetation_dataset)
}

get_clean_tile <- function(brick){
  #spectra = readr::read_csv("/Users/sergiomarconi/Documents/GitHub/hiPyRneon/outdir/spectra/itc2px_reflectance.csv")
  spectra = spectra %>% dplyr::select(matches("individualID|band"))
  foo=spectra[-1]
  mask = foo > 10000
  foo[mask] <-  NA
  spectra[-1]=foo
  spectra = dplyr::na_if(spectra, -9999)
  spectra = spectra[complete.cases(spectra),]
  spectra=optical.filter(spectra)
  spectra = spectra[complete.cases(spectra),]
  spectra = spectra[c(1,15:365)]
  spectra[-1] = hiper_features(spectra[-1])
  #spectra[-1] = scale(spectra[-1])
  #vst = sf::read_sf("./indir/traits_polygon_dataset.shp") %>%
  # dplyr::select(indvdID, domanID, siteID, taxonID, scntfcN, elevatn, nlcdCls, cllctDt, ntrgnPr, 
  #               crbnPrc, CNratio , extrcCC, extrCAC, extrCBC, dryMass, dryMssF, 
  #               lfMssPA, leafAre, lgnnPrc, clllsPr, plntStt)
  
  vst = sf::read_sf("./indir/NEON_polygon_dataset.shp") %>%
    dplyr::select(indvdID, domanID, siteID, taxonID, scntfcN)  %>% filter(siteID == "HARV")
  colnames(vst)[1] = "individualID"
  #vst$cllctDt = lubridate::month(vst$cllctDt)
  vegetation_dataset = inner_join(vst, spectra) %>% unique %>% data.frame
  #vegetation_dataset = vegetation_dataset[complete.cases(vegetation_dataset),]
  #readr::write_csv(vegetation_dataset, "~/Documents/Data/vegetation_dataset.csv")
  return(vegetation_dataset)
}

