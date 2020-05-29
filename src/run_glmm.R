
run_glmm <- function(token){
  library(reticulate)
  library(tidyverse)
  library(brms)
  
  dataset = readr::read_csv("./outdir/brdf_corrected_kld_hist.csv")
  cfc_reduced_dims = readRDS("//orange/ewhite/s.marconi/Chapter3/hdr_cfc.rds")
  #cfc_data = cbind(dataset[c(1:21, 374:376, 378:380)], cfc_reduced_dims$firstpc) #, 374:380
  cfc_data = dataset
  IDs = readr::read_csv("./indir/cfc_higher_3.csv")
  cfc_data = cfc_data %>% filter(individualID %in% unlist(IDs))
  nComp = cfc_reduced_dims$bands_grouping %>% unique %>% length
  refl_fixed_eff = paste("(", colnames(cfc_data)[-c(1:4, 6:15)],")", collapse = " + ", sep = "")
  tr_nm = colnames(dataset)[c(6:15)]     #c(9,10,12, 13,14,17, 19, 20)]
  #tr_nm=tr_nm[1]
  f1 <- paste(paste("mvbind(", paste(tr_nm, collapse = " , "), ") ~ ", sep = ""),
              refl_fixed_eff," + ", 
              paste("(1 | siteID)"))
  
  #cfc_data[9:20] = log(cfc_data[9:20])
  train = cfc_data %>% group_by(individualID) %>%
    sample_n(1)
  set.seed(1987)
  train = train %>% group_by(siteID) %>% sample_frac(0.8)
  train = cfc_data %>% filter(individualID %in% unique(train$individualID))
  test = cfc_data %>% filter(!individualID %in% unique(train$individualID))
  train["elevation"] = train["elevation"] %>%
    mutate_if(is.numeric, scale)
  #train[,c(6, 23:45)] = train[,c(6, 23:45)] %>%
  #  mutate_if(is.numeric, scale)
  
  scaling_fct = list()
  for(tr in colnames(train["elevation"])){
    scaling_fct[[tr]] = cbind.data.frame(tr, attr(train[[tr]], "scaled:center"), 
                                attr(train[[tr]], "scaled:scale"))
  }
  scaling_fct = do.call(rbind.data.frame, scaling_fct)
  colnames(scaling_fct) <- c("feature", "center", "scale")
  
  test["elevation"] = lapply(c(1:nrow(scaling_fct)), function(x) {
    dat = test["elevation"]
    scale(dat[,x], center = scaling_fct[x,"center"], scale = scaling_fct[x,"scale"])})
  set.seed(token)
  train = train %>% group_by(individualID) %>%
    sample_n(1)
  
  #train[,c(9:20)] = log(train[,c(9:20)])
  fit <- brm(f1,# + set_rescor(FALSE),  
             data = train,
             cores = 2,
             seed = 12,
             family = lognormal(),
             control = list(adapt_delta = 0.99),
             thin = 10, #refresh = 0,
             prior = set_prior(horseshoe(df = 3)),
             chains = 2,
             iter = 3000)
  
  #test[,c(9:20)] = log(test[,c(9:20)])
  #print(bayes_R2(fit))
  test_r2 = bayes_R2(fit, newdata = test)
  #
  prds = predict(fit, newdata = test)
  
  #prds = prds[,1,]
  derived = lapply(1:4, function(x){
    quart = prds[,x,] %>% data.frame
    quart["individualID"] = test$individualID
    quart  %>% group_by(individualID) %>%
      summarize_all(median)
  })
  get_mod_r2 <- function(pred, obs){
    1 - sum((pred - obs)^2) / sum((obs - mean(obs, na.rm=T))^2)
  } 
  
  
  test_itc = test %>% select(c("individualID", tr_nm)) %>% unique
  
  R2 = lapply(1:length(tr_nm), function(x){
    data_test = inner_join(derived[[1]][,c(1, x+1)], 
                           test_itc[,c(1, x+1)], by="individualID")
    get_mod_r2((data_test[[2]]), (data_test[[3]]))
  })
  R2 = unlist(R2)
  R2
  #valid = kfold(fit, K = 4, cores = 2, folds = "stratified", group = "siteID")
  
  res = list(mod = fit, test_set = test, itcR2 = R2, br2 = test_r2, scaling = scaling_fct)
  saveRDS(res, paste("./outdir/mods/", token, "_md.rds", sep=""))
  return(R2)
}

library(parallel)
cl <- makeCluster(22)
stack_r2 = parLapply(cl, 1:300, run_glmm)
stopCluster(cl)
