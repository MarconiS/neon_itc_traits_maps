#read shp and df and append them
library(sf)
library(tidyverse)

shp = read_sf("/Users/sergiomarconi/Documents/GitHub/Chapter3/traits_bboxes/cfc_polygond.shp")
colnames(shp)[1:2] <- c("uncFlag", "individualID")
df = read_csv("/Users/sergiomarconi/Documents/GitHub/neonVegWrangleR/outdir/field_data/cfc_with_height.csv")
df =  df %>% filter(height_smpl > 3) 
df$individualID = substring(df$individualID, 14)  
shp = inner_join(shp, df, by="individualID")
colnames(shp)
shp = shp %>% select(uncFlag, individualID, eventID, domainID, siteID, plotID, subplotID, 
                     taxonID, scientificName, utmZone, elevation, nlcdClass, shpEasting, shpNorthing, plotLatitude, plotLongitude, 
                     longitude, latitude, collectDate, nitrogenPercent, carbonPercent, CNratio, 
                     extractCarotConc, extractChlAConc, extractChlBConc, dryMass, dryMassFraction, 
                     leafMassPerArea, leafArea, ligninPercent, cellulosePercent, plantStatus)
shp = unique(shp)

sf::write_sf(shp, "./indir/traits_polygon_data_may.shp")
shp$individualID %>% table %>% sort %>% tail(10)




#read shp and df for harv
library(sf)
library(tidyverse)

#shp = read_sf("/Users/sergiomarconi/Documents/Data/annotations/HARV.shp")
shp  = sf::read_sf("./weak_label/indir/shp/HARV.shp")
shp$treeID = paste("HARV", shp$treeID, sep=".")
colnames(shp)[1:2] <- c("uncFlag", "individualID")
df = read_csv("./weak_label/indir/csv/vst_field_data.csv")

df$individualID = substring(df$individualID, 14)  
shp = inner_join(shp, df, by="individualID")
colnames(shp)
shp = shp %>% select(uncFlag, individualID, eventID, siteID, plotID, 
                     taxonID, itcLongitude, itcLatitude, plantStatus)
shp = unique(shp)
shp = shp %>% group_by(individualID)%>% sample_n(1)
sf::write_sf(shp, "./weak_label/indir/shp/HARV_full_preliminary.shp")
shp$individualID %>% table %>% sort %>% tail(10)

df$canopyPosition %>% unique




