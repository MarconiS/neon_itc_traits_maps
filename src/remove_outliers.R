remove_outliers <- function(spectra, method = "lof"){
  library(bigutilsr)
  pca <- prcomp(spectra[-1], scale. = TRUE)
  U <- pca$x
  # library(ggplot2)
  # theme_set(bigstatsr::theme_bigstatsr(0.8))
  # qplot(U[, 1], U[, 2]) + coord_equal()
  # 
  if(method == "pca_dist"){
    ids <- apply(U, 2, function(x) which( (abs(x - median(x)) / mad(x)) > 6 )) %>%
      Reduce(union, .)
  }
  if(method == "lof"){
    llof <- LOF(U)
    div_lim = hist_out(llof)$lim[2]
    ids = which(llof > div_lim)
  }

  eigval <- pca$sdev^2
  npc = pca_nspike(eigval)
  return(list(outliers = ids, npc = npc))
}
