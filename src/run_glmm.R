
run_glmm <- function(token){
  library(reticulate)
  library(tidyverse)
  library(brms)
  
  dataset = readr::read_csv("./indir/brdf_18_cfc_june.csv")
  dataset = dataset[complete.cases(dataset),]
  dataset = dataset %>% dplyr::select(-one_of("taxonID.y","siteID.y","domainID.y" ))
  colnames(dataset)[2:4] = c("taxonID","siteID","domainID" )
  cfc_data = dataset
  IDs = readr::read_csv("./indir/brdf_18_cfc_june.csv")
  IDs = dataset %>% filter(!siteID %in%c("SERC", "ORNL")) %>%
    #filter(height >2) %>%
    #filter(stemDiameter >5) %>% 
    select(individualID) %>% unique

  cfc_data = cfc_data %>% filter(individualID %in% unlist(IDs))
  refl_fixed_eff = paste("(", colnames(cfc_data)[-c(1:5, 8:17)],")", collapse = " + ", sep = "")
  tr_nm = colnames(dataset)[c(8:17)]     #c(9,10,12, 13,14,17, 19, 20)]
  f1 <- paste(paste("mvbind(", paste(tr_nm, collapse = " , "), ") ~ ", sep = ""),
              refl_fixed_eff," + ", 
              "(1 | dm | domainID) + ",
              paste("(1 | gr | siteID)"))
  
  train = cfc_data %>% group_by(individualID) %>%
    slice(1)

  set.seed(1987)
  train = train %>% group_by(siteID) %>% sample_frac(0.7)
  train = cfc_data %>% filter(individualID %in% unique(train$individualID))
  test = cfc_data %>% filter(!individualID %in% unique(train$individualID))
  set.seed(1987)
  oob = test %>% select(individualID, siteID) %>% unique %>% group_by(siteID) %>% sample_frac(0.4)
  oob = test %>% filter(individualID %in% unique(oob$individualID))
  test = test %>% filter(!individualID %in% unique(oob$individualID))
  
  train[-c(1:5, 8:17)] = train[-c(1:5, 8:17)] %>%
    mutate_if(is.numeric, scale)

  scaling_fct = list()
  for(tr in colnames(train[-c(1:5, 8:17)])){
    scaling_fct[[tr]] = cbind.data.frame(tr, attr(train[[tr]], "scaled:center"), 
                                         attr(train[[tr]], "scaled:scale"))
  }
  scaling_fct = do.call(rbind.data.frame, scaling_fct)
  colnames(scaling_fct) <- c("feature", "center", "scale")
  
  test[-c(1:5, 8:17)] = lapply(c(1:nrow(scaling_fct)), function(x) {
    dat = test[-c(1:5, 8:17)]
    scale(dat[,x], center = scaling_fct[x,"center"], 
          scale = scaling_fct[x,"scale"])})
  oob[-c(1:5, 8:17)] = lapply(c(1:nrow(scaling_fct)), function(x) {
    dat = oob[-c(1:5, 8:17)]
    scale(dat[,x], center = scaling_fct[x,"center"], 
          scale = scaling_fct[x,"scale"])})
  set.seed(token)
  train = train %>% group_by(individualID) %>%
    sample_n(1)
  #set priors
  ## check which parameters can have priors
  #get_prior(f1, data = train, family = lognormal())
  
  ## define some priors          
  prior <- c(set_prior("horseshoe(3)", class = "b"))
             #, set_prior("cauchy(0,2)", class = "sd" 
             #,  group = "subject", coef = "Intercept"))
  
  
  fit <- brm(f1,
             data = train,
             cores = 2,
             seed = 12,
             family = lognormal(),
             control = list(adapt_delta = 0.99),
             thin = 10, 
             prior = set_prior(horseshoe()),
             chains = 2,
             iter = 4000)
    #kfold <- kfold(fit, chains = 1, folds = "stratified", group ="siteID", K = 3, save_fits = T)
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
    
    res = list(mod = fit, test_set = test, itcR2 = R2, br2 = test_r2, scaling = scaling_fct, oob_data = oob)
    saveRDS(res, paste("./outdir/mods/", token, "_md.rds", sep=""))
  return(R2)
}

library(parallel)
cl <- makeCluster(2)
stack_r2 = parLapply(cl, sample.int(10000, 2), run_glmm)
stopCluster(cl)
