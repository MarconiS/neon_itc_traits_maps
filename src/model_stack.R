ls_mod = list.files("./outdir/mods/", full.names = T)
R2s <- list()
for(pt in 1:length(ls_mod)){
  tryCatch({
    mod = readRDS(ls_mod[pt])
    R2s[[pt]] = mod$itcR2
  })
}
histsr2 = do.call(rbind.data.frame, R2s)
colnames(histsr2) <- c("N", "C",  "CC", "CAC", "CBC",  "dryMassFrac", "LMA", "lign", 'cell')
histsr2[["sumR2"]] = apply(histsr2[,c(1,2,7,8)], 1, sum)
foo = histsr2[order(histsr2$sumR2),]
foo = tail(foo, 250)

mods_to_save = rownames(foo)
mods = list()
for(pt in as.numeric(mods_to_save)){
  ft = readRDS(ls_mod[pt])
  mod[[pt]] = ft$mod
}
library(brms)
md_weights = model_weights(mod[[1]],mod[[5]],mod[[9]],mod[[13]],mod[[17]],
                           mod[[2]],mod[[6]],mod[[10]],mod[[14]],mod[[18]],
                           mod[[3]],mod[[7]],mod[[11]],mod[[15]],mod[[19]],
                           mod[[4]],mod[[8]],mod[[12]],mod[[16]],mod[[20]], 
                           mod[[21]],mod[[22]],mod[[23]],mod[[24]],mod[[25]],
                           mod[[26]],mod[[27]],mod[[28]],mod[[29]],mod[[30]],
                           mod[[31]],mod[[32]],mod[[33]],mod[[34]],mod[[35]],
                           mod[[36]],mod[[37]],mod[[38]],mod[[39]],mod[[40]],
                           mod[[41]],mod[[42]],mod[[43]],mod[[44]],mod[[45]],
                           mod[[46]],mod[[47]],mod[[48]],mod[[49]],mod[[50]])

#ideally check all and get only models whose weight sum to 99+
md_weights = md_weights %>% sort(decreasing = T)
thr = 0
md_choose = NULL
for(l in 1:length(md_weights)){
  thr = thr  + md_weights[l]
  md_choose = c(md_choose, names(md_weights[l]))
  if(thr > 0.99){
    break
  }
}

#md_choose = substring(md_choose, 6)  
#md_choose = do.call(c, strsplit(md_choose, split = "]]"))
test_data = ft$test_set
mu = sd = rep(0, nrow(test_data))
for(ii in 1:length(md_choose)){
  tmp = predict(mod[[as.integer(md_choose[ii])]], newdata = test_data,
                nsamples =200, robust = T)
  #expected value adn model error
  mu = mu + tmp[,1,] * md_weights[ii]
  sd = sd + tmp[,2,] * md_weights[ii]
}
mu = data.frame(mu)
colnames(mu) = c("N", "C",  "CN", "CC", "CAC", "CBC", "LMA", "lign", 'cell')
mu['individualID'] = test_data$individualID
mu = mu %>% group_by(individualID) %>%
  summarize_all(median)

yobs = test_data %>% select(individualID, ntrgnPr, crbnPrc, CNratio, 
                            extrcCC, extrCAC, extrCBC,
                            lfMssPA, lgnnPrc, clllsPr) %>%
  group_by(individualID) %>% summarize_all(median)

#calculate R2
R2 = list()
for(ii in 2:ncol(mu)){
  R2[[ii]] = get_mod_r2(mu[[ii]], yobs[[ii]])
}
unlist(R2)