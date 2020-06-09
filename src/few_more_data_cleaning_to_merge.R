
#clean data
#vegetation_dataset <- read_csv("ccbid/support_files/vegetation_dataset.csv")
vegetation_dataset = get_full_dataset(final_dataset)
vegetation_dataset = vegetation_dataset %>% filter(band_312 <0.025)  
cnd1 = (vegetation_dataset[24:60] > 0.04)
idx <- which(apply(cnd1, 1, any))
if(length(idx) !=0){
  vegetation_dataset = vegetation_dataset[-idx,]
}
cnd1 = (vegetation_dataset[100:200] > 0.125)
idx <- which(apply(cnd1, 1, any))
if(length(idx) !=0){
  vegetation_dataset = vegetation_dataset[-idx,]
}

outlr <- vegetation_dataset %>%
  dplyr::select(matches("band|individualID")) %>%
  remove_outliers()
#plot_spectra(vegetation_dataset[-outlr$outliers,])

write_csv(vegetation_dataset[-outlr$outliers,], "./outdir/harv_reflectance.csv")


# list_ids = df %>% filter(canopyPosition %in% c("Full sun",  "Partially shaded", "Open grown")) %>%
#   filter(siteID == "HARV")
# list_ids =  substring(list_ids$individualID, 19)  
#vegetation_dataset = vegetation_dataset %>% dplyr::filter(individualID %in% list_ids)
#get species classification format dataset
sc_names = lapply(1:nrow(vegetation_dataset), function(x) 
  strsplit(vegetation_dataset[["scntfcN"]][x], split = " ")[[1]][1:2])
sc_names = do.call(rbind.data.frame, sc_names)
colnames(sc_names) <- c("species", "genus")
vegetation_dataset <- vegetation_dataset %>% dplyr::select(-one_of("scntfcN"))
vegetation_dataset <- cbind.data.frame(vegetation_dataset, sc_names)
genus_id =  substr(vegetation_dataset[["taxonID"]], start = 1, stop = 2) 
vegetation_dataset <- cbind.data.frame(vegetation_dataset, genus_id)

vegetation_dataset[["individualID"]] = factor(vegetation_dataset[["individualID"]])
too_rare = vegetation_dataset %>% group_by(individualID) %>% top_n(1,wt = band_14) %>% 
  ungroup %>% dplyr::select(taxonID, siteID) %>% table %>% data.frame


# #remove what is not a species, and correct typos
remove = c("2PLANT", "PINUS", "OSVI", "CAGL8", "FAGR","QUVE")

vegetation_dataset[vegetation_dataset$taxonID =="BEPAP", "taxonID"] = "BEPA"
vegetation_dataset[vegetation_dataset$taxonID =="ACSAS", "taxonID"] = "ACSA3"
vegetation_dataset[vegetation_dataset$taxonID =="PICEA", "taxonID"] = "PIAB"

# change too few in OTHR
for(ii in unique(too_rare$siteID)){
  which_rare = too_rare %>% filter(siteID == ii) %>% filter(Freq <6) %>% filter(Freq > 0)
  vegetation_dataset[vegetation_dataset$taxonID %in% which_rare$taxonID &
                       vegetation_dataset$siteID %in% which_rare$siteID, "taxonID"] = "OTHR"
}
vegetation_dataset[vegetation_dataset$taxonID %in% remove, "taxonID"] = "OTHR"



#remove other
vegetation_dataset = vegetation_dataset %>% filter(taxonID != "OTHR")
#split data
test_itcs = vegetation_dataset %>% group_by(siteID, taxonID) %>% 
  dplyr::select(individualID) %>% unique %>%
  dplyr::sample_frac(size = 0.3)
train_set = vegetation_dataset %>% filter(!(individualID %in% test_itcs$individualID))
test_set= vegetation_dataset %>% filter(individualID %in% test_itcs$individualID)

#mutate crown ids into number for numpy compatibility
test_set[["individualID"]] =as.numeric(test_set[["individualID"]])
train_set[["individualID"]] =as.numeric(train_set[["individualID"]])

vst_species_train = train_set %>% 
  dplyr::select(individualID, taxonID, genus_id, species, genus, siteID)
vst_species_test = test_set %>% 
  dplyr::select(individualID, taxonID, genus_id, species, genus, siteID)
train_set = train_set %>% 
  dplyr::select(matches("band|individualID"))
test_set = test_set %>% ungroup %>% 
  dplyr::select(matches("band|individualID"))

colnames(vst_species_test)[1:2] = colnames(vst_species_train)[1:2] = c("crown_id", "species_id")
colnames(train_set)[1] = colnames(test_set)[1] = c("crown_id")

#save dataset to be used in stanford ccb
readr::write_csv(train_set, "./outdir/cfc_sp_training.csv")
readr::write_csv(test_set,"./outdir/cfc_sp_testing.csv")
readr::write_csv(unique(vst_species_train), "./outdir/cfc_sp_species_id.csv")
readr::write_csv(unique(vst_species_test), "./outdir/cfc_sp_test_id.csv")

vst_species_train %>% unique %>% dplyr::select(siteID) %>% table

