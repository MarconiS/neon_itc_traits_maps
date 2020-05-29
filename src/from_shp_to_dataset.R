library(tidyverse)
library(sf)
#read shp file
from_shp_to_ml_dataset <- function(shp = NULL, files_path = NULL, outpath = "./outdir/HARV_geo_reflectance.csv", bff = 2, save_imgs = F){
  shp = sf::read_sf("./indir/traits_polygon_data_may.shp")
  colnames(shp)[2] <- "individualID"
  # shp  =read_sf("../Chapter1/corrected_res/OSBS_crown_polygons/OSBS_sample_polygons_edits_Feb2018.shp")
  #shp  = sf::read_sf("./weak_label/indir/shp/HARV_full_preliminary.shp") #"~/Documents/Data/gfor_HARV.shp") #"./weak_label/indir/shp/HARV_full_preliminary.shp")
  #shp = sf::read_sf("/Users/sergiomarconi/Documents/Data/NEON/VST/vst_top_canopy.csv")
  #shp = sf::st_as_sf(shp, coords = c("itcLongitude", "itcLatitude"), crs = 4326)
  files_path =  "~/Documents/Data/RS/plots/bdrf/" #"../neonVegWrangleR/outdir/plots/hsi/" #"~/Documents/Data/RS/plots/hsi/"
  list_hsi = list.files(files_path, pattern = ".tif", full.names = T)
  
  #loop thrpugh each id and get both clip and append 
  df = list()
  for(ii in 1:nrow(shp)){
    itc = shp[ii,]
    fls = grep(shp$plotID[ii], list_hsi, value=TRUE)
    tmp = list()
    for(pltID in fls){
      tryCatch({
        plt = raster::brick(pltID)
        itc = st_transform(itc, crs = raster::crs(plt))
        #itc = st_centroid(itc)
        if(st_geometry_type(itc) == "POLYGON"){
          dat = raster::extract(plt, itc) %>% data.frame
        }else{
          dat = raster::extract(plt, itc, buffer = 1) %>% data.frame
        }
        if(length(dat) !=0 & any(!is.na(dat))){
          if(ncol(dat) < 369){
            dat = t(dat)
          }
          colnames(dat) = paste("band", 1:369, sep="_")
          tmp[[pltID]] = cbind.data.frame(shp$individualID[ii], dat)
          if(save_imgs == T){
            #   #extract and save a clipped image
            box = raster::extent(itc)
            box[c(1,3)] = box[c(1,3)]-bff
            box[c(2,4)] = box[c(2,4)]+bff
            crp = raster::crop(plt, box)
            # # raster::plotRGB(crp, 18, 55, 87, stretch="lin")
            raster::writeRaster(crp, paste("./outdir/crops/brdf/", shp$individualID[ii], 
                                           ".tif", sep=""),overwrite=TRUE)
            
            rgbID =gsub("\\<hsi\\>","mage",pltID)
            rgb  = raster::brick(rgbID)
            crp = raster::crop(rgb, box)
            # # raster::plotRGB(crp, 1, 2, 3, stretch="lin")
            raster::writeRaster(crp, paste("./outdir/crops/rgb/", shp$individualID[ii], 
                                           ".tif", sep=""),overwrite=TRUE)
            
            #
            if(length(tmp) == 0){
              warning(paste(shp$plotID[ii], "missing"))
            }else{
              rownames(tmp[[pltID]]) = NULL
              colnames(tmp[[pltID]])[1] = "individualID"
            }
          }
        }
      })
    }
    tmp = do.call(rbind.data.frame, tmp)
    df[[ii]] = tmp
  }
  final_dataset = do.call(rbind, df)
  write_csv(final_dataset, outpath)
  return(final_dataset)
}
