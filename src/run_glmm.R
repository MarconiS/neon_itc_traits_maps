
run_glmm <- function(token){
  library(reticulate)
  library(tidyverse)
  library(brms)
  
  dataset = readr::read_csv("./indir/brdf_june.csv")
  dataset = dataset[complete.cases(dataset),]
  cfc_data = dataset
  IDs = readr::read_csv("./indir/cfc_height_andstem.csv")
  IDs = IDs %>% filter(!siteID %in%c("SERC", "ORNL")) %>%select(individualID) %>% unique

  cfc_data = cfc_data %>% filter(individualID %in% unlist(IDs))
  refl_fixed_eff = paste("(", colnames(cfc_data)[-c(2:15)],")", collapse = " + ", sep = "")
  tr_nm = colnames(dataset)[c(6:15)]     #c(9,10,12, 13,14,17, 19, 20)]
  f1 <- paste(paste("mvbind(", paste(tr_nm, collapse = " , "), ") ~ ", sep = ""),
              refl_fixed_eff," + ", 
              paste("(1 | gr | siteID)"))
  
  train = cfc_data %>% group_by(individualID) %>%
    slice(1)

  set.seed(1987)
  train = train %>% group_by(siteID) %>% sample_frac(0.8)
  train = cfc_data %>% filter(individualID %in% unique(train$individualID))
  test = cfc_data %>% filter(!individualID %in% unique(train$individualID))
  train["elevation"] = train["elevation"] %>%
    mutate_if(is.numeric, scale)

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
  
  fit <- brm(f1,
             data = train,
             cores = 2,
             seed = 12,
             family = lognormal(),
             control = list(adapt_delta = 0.99),
             thin = 10, 
             #prior = set_prior(horseshoe()),
             chains = 2,
             iter = 3000)
  
  test_r2 = bayes_R2(fit, newdata = test)
  #
  prds = predict(fit, newdata = test)
  
  res = list(mod = fit, test_set = test, itcR2 = R2, br2 = test_r2, scaling = scaling_fct)
  saveRDS(res, paste("./outdir/mods/", token, "_md.rds", sep=""))
  return(R2)
}

library(parallel)
cl <- makeCluster(32)
stack_r2 = parLapply(cl, sample.int(10000, 40), run_glmm)
stopCluster(cl)
