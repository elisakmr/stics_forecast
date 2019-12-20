load(file = file.path(dir_data, "extract", "years.RData"))

# generate random year series
 for (b in 1:20){
    rand_year[[b]] <- rand_year[[b]]-1992
  }
  save(rand_year, file = file.path(dir_data, "extract", "years.RData"))

  tercquint("obs","seasonal",array_obs,24,34,24)[1,,3]
  seasbin_ref[1,,3]
  seas_obs<-apply(array_obs[1,,],1, sum)
  qsup <- quantile(as.vector(seas_obs), 2/3)
  which(seas_obs>qsup)



  veriApply("EnsRoca", fcst=apply(ps_original,c(1,2,3),sum), obs=apply(ref_original,c(1,2),sum), prob=1:2/3, parallel = TRUE)$cat1
  EnsRoca(apply(ps_original[18,,,],c(1,2),sum), apply(ref_original[18,,],c(1),sum))

  veriApply("EnsRoca", fcst=apply(ps_original[18,,,],c(1,2),sum), obs=apply(ref_original[18,,],c(1),sum), prob=1:2/3, parallel = TRUE)$cat3



indicateur <- "DRAINIRRIG"

load(file = file.path(dir_data, "extract", "years.RData")) # list of years in each bootstrap iteration

nom_ref <- paste0(indicateur,1,"_ref.RData")
load(file = file.path(dir_data, "extract", nom_ref))
ngrids <- dim(array_obs)[1]

ps_original <- array(NA, dim=c(ngrids, nan, nrun, ndimT))
ref_original <- array(NA, dim=c(ngrids, nan, ndimT))

for (year in an){
  y <- which(an==year)
  # bootstrap iteration index of the year selected

  list_match <- lapply(rand_year, function(x) which(year == x)) # list of matchs between list of years randomly generated for the bootstrap and the year we are interested in

  y_index <- list_match[[1]] # year index within the bootstrap iteration
  b_index <- 1 # bootstrap iteration index
  for (b in 1:(nboot-1)){ # pick the first iteration and first index matching the year we are looking for
    if (length(list_match[[b]])<1 & length(y_index)<1) {y_index <- list_match[[b+1]]; b_index <- b+1}
    else {y_index <- y_index; b_index <- b_index}
  }

  nom_prev <- paste0(indicateur,b_index,"_prev.RData")
  load(file = file.path(dir_data, "extract", nom_prev))
  ps_original[,y,,] <- array_prev[,y_index[1],,]

  nom_ref <- paste0(indicateur,b_index,"_ref.RData")
  load(file = file.path(dir_data, "extract", nom_ref))
  ref_original[,y,] <- array_obs[,y_index[1],]
}

apply(ref_original[18,,],1,sum)


