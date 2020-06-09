get_climate <- function(field_data = NULL, tmppath = './tmp/'
                        , startdate = 1998
                        , enddate = 2018
                        , provider = "daymet"
){
  if(provider=="daymet"){
    download_point_daymet(field_data = field_data,
                          outpath = tmppath,
                          startdate=startdate,
                          enddate = enddate)
    melt_daymet_trends(path = tmppath)
    climate  = readr::read_csv("./indir/climate_features.csv")
    return(climate)
  }
}

#' from coordinates, extract data for each of the data-products of interest
#'
#'
#' @return
#' @export
#' @import daymetr lubridate
#' @examples
#' @importFrom magrittr "%>%"

melt_daymet_trends <- function(path = "./tmp/"
                               , outfile = "./indir/climate_features.csv"
){
  
  dataset = data.frame(matrix(NA, ncol = 9, nrow = 0))
  colnames(dataset) <- c("month","ts_daylength", "ts_prec","ts_rad",
                         "ts_melt","ts_tmax","ts_tmin","ts_vp", "individualID")
  ls_dat = list.files(path, pattern = "csv")
  for(ii in ls_dat[-1]){
    id_clim <- read.csv(paste(path,ii, sep="/"), skip=7)
    colnames(id_clim) <- c("year", "month", "daylength", "prec", "srad", "snow_melt",
                           "tmax", "tmin", "vp")
    id_clim$month <- as.Date(id_clim$month-1, origin = "1995-01-01")# %>%
  #@month(id_clim$month)
    point_features <- id_clim %>%
      dplyr::group_by(month, year) %>%
      dplyr::summarise(daylength = mean(daylength), prec = sum(prec),
                       srad = mean(srad), snow_melt = sum(snow_melt),
                       tmax = max(tmax), tmin = min(tmin), vp = mean(vp))
    
    point_features <- point_features %>%
      dplyr::group_by(month) %>%
      dplyr::summarise(daylength = mean(daylength), prec = mean(prec),
                       rad = mean(srad), snow_melt = mean(snow_melt),
                       tmax = mean(tmax), tmin = mean(tmin), vp = mean(vp))
    point_features[["individualID"]] <- gsub('.{14}$', '', ii)
    dataset = rbind(dataset, point_features)
  }
  #write_csv(dataset, "./DMT_retriever/climate_features.csv")
  dataset$month =  lubridate::month(dataset$month)
  library(data.table) ## v >= 1.9.6
  clim_dat <- dcast(melt(as.data.table(dataset), id.vars=c("individualID", "month")),
                    individualID~variable+month, fun.aggregate = mean)
  readr::write_csv(clim_dat, outfile)
  return(clim_dat)
}

#' from coordinates, extract data for each of the data-products of interest
#'
#'
#' @return
#' @export
#' @examples
#' @import tidyverse
#' @importFrom magrittr "%>%"
download_point_daymet <- function(listSites = NULL,
                                  field_data = NULL,
                                  outpath = './tmp/',
                                  startdate=1998,
                                  enddate= 2018){
  #field_data <- dplyr::select(field_data, individualID, latitude, longitude)
  if(is.null(field_data[["latitude"]]) || is.na(field_data[["longitude"]])){
    new_dat <- get_lat_long(field_data)
    daymet_coords <- cbind(as.character(field_data[["individualID"]]),
                           new_dat[["northing"]],
                           new_dat[["easting"]]) %>% unique
  }else{
    daymet_coords <- cbind(as.character(field_data[["individualID"]])
                           , field_data[["latitude"]]
                           , field_data[["longitude"]]) %>%
      unique
  }
  readr::write_csv(data.frame(daymet_coords), './tmp/daymet_metadata.csv')
  
  library(daymetr)
  download_daymet_batch(file_location = './tmp/daymet_metadata.csv',
                        start = startdate,
                        end = enddate,
                        internal = F,
                        path = outpath)
}

predict_data_pca_transforamtion <- function(climate_trends, pca_mod) {
  #PCA by variable
  var <- c("daylength", "prec", "rad", "snow_melt", "tmax", "tmin", "vp")
  climate_pca = NULL
  climate_features = NULL
  # climate_trends[-1] = scale(climate_trends[-1], center = scaling_vars$center[-c(1:7)],
  #                           scale = scaling_vars$scale[-c(1:7)])
  
  for(ii in var){
    climate_features[[ii]] <- climate_trends %>%
      select(contains(ii)) 
    pca_mod[[ii]]$var$coord = data.frame((pca_mod[[ii]]$var$coord))
    climate_pca[[ii]]  <-  FactoMineR::predict.PCA(pca_mod[[ii]], climate_features[[ii]])
    #get pca transformed data
    climate_features[[ii]] = climate_pca[[ii]]["cos2"]
  }
  climate_features = do.call(cbind.data.frame, climate_features)
  colnames(climate_features) <- var
  climate_features <- data.frame(climate_trends[1], climate_features)
  colnames(climate_features)[1] <- "individualID"
  return(list(climate_pca=climate_pca, climate_features=climate_features))
}


field_data = df %>% select(individualID, longitude, latitude) %>% unique
climate_FIA <- readRDS("../Chapter2/indir/climate_FIA.rds")
#climate_trends = get_climate(field_data)
climate_trends = readr::read_csv("./indir/climate_features.csv")
climate <- predict_data_pca_transforamtion(climate_trends, 
                                           climate_FIA$climate_pca)
climate_features = climate$climate_features


