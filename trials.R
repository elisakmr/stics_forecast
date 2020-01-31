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

#################################################################

                                               ## GRIDDED DET SCORES ##


acc_prev <- vector()
msess <- vector()

nom_ref <- paste0(indicateur,b,"_ref.RData")
nom_prev <- paste0(indicateur,b,"_prev.RData")
nom_clim <- paste0(indicateur,b,"_clim.RData")
nom_an_prev <- paste0(indicateur,b,"_anomaly_prev.RData")
nom_an_clim <- paste0(indicateur,b,"_anomaly_clim.RData")
nom_an_obs <- paste0(indicateur,b,"_anomaly_obs.RData")

load(file = file.path(dir_data, "extract", nom_ref))
load(file = file.path(dir_data, "extract", nom_prev))
load(file = file.path(dir_data, "extract", nom_clim))
load(file = file.path(dir_data, "extract", nom_an_prev))
load(file = file.path(dir_data, "extract", nom_an_clim))
load(file = file.path(dir_data, "extract", nom_an_obs))

ngrids <- dim(array_obs)[1]
#ndimT <- dim(array_obs)[3]

# ensemble mean
prev_mean <- apply(array_prev[,,,], c(1,2,4), mean)
clim_mean <- apply(array_clim[,,,], c(1,2,4), mean)


### SEASONAL ###


# BE CAREFUL ON AGREGATION METHOD!!!!
det_seas <- array(NA, dim = c(ngrids,3))

prev_seas <- apply(array_prev[,,,], c(1,2,3), sum) #### aggregation method to be changed!!
clim_seas <- apply(array_clim[,,,], c(1,2,3), sum)
ref_seas <- apply(array_obs[,,], c(1,2), sum)

prev_seas_mean <- apply(prev_seas, c(1,2), mean)
clim_seas_mean <- apply(clim_seas, c(1,2), mean)

anoprev_seas <- apply(prev_anomaly, c(1,2,3), sum)
anoref_seas <- apply(obs_anomaly, c(1,2), sum)

# MSESS

msess <- veriApply("EnsMsess", fcst=t(prev_seas_mean), fcst.ref=t(clim_seas_mean), obs=t(ref_seas), parallel = TRUE)

# ACC #

acc_prev <- veriApply("EnsCorr", fcst=anoprev_seas, obs=anoref_seas)
det_seas <- array(c(msess, acc_prev), dim=c(ngrids,2))

nom_prevseas <- paste0(indicateur, b,"scoredetseas_grid_prev.RData")

save(det_seas, file = file.path(dir_data, "score", nom_prevseas))


                                              ## BOX DET SCORES ##


coeff <- vector()
seas_ps <- matrix(ncol = 2, nrow = 6)

  load(file = file.path(dir_data, "extract", paste0("lat_",indicateur,".RData")))
  boot_seas <- list()

    # loading gridded scores

    # seasonal
    ngrids <- dim(det_seas)[1]

    val_prevseas <- list()

    for (i in 1:ngrids){
      val_prevseas[[i]] <- det_seas[i,]*cos(pi*lat[i]/180)
      coeff[i] <- cos(pi*lat[i]/180)
      # removing NAs
      if (length(which(is.na(det_seas[i,])))>0){
        coeff[i] <- NA
        val_prevseas[[i]] <- rep(NA,2)
      }
    }

    boot_seas[[b]] <- apply(array(unlist(val_prevseas) , c(2,ngrids)), 1, sum, na.rm = T)/sum(coeff, na.rm = T)




                                               ### PROB SCORE EXPLORING ###
## obs

    nom_ref <- paste0(indicateur,b,"_ref.RData")
    nom_prev <- paste0(indicateur,b,"_prev.RData")

    load(file = file.path(dir_data, "extract", nom_ref))
    load(file = file.path(dir_data, "extract", nom_prev))

    seas_obs <- apply(array_obs[,,9:20],c(1,2), sum)
    seas_prev <- apply(array_prev[,,,9:20],c(1,2,3), sum)
    ps_season <- apply(seas_prev, 2, mean)
    ref_season <- apply(seas_obs, 2, mean)

    # ps_season_int <- apply(ps_original[,,,9:20], c(1,2,3), sum)
    # ps_season <- apply(ps_season_int, 2, mean)
    # ref_season_int <- apply(ref_original[,,9:20], c(1,2), sum)
    # ref_season <- apply(ref_season_int, 2, mean)

    df_raw <- data.frame()
    for (y in 1:nan){

      df_raw <- rbind(df_raw, data.frame(type = "prevision", year=an[y], value = ps_season[y]))
      df_raw <- rbind(df_raw, data.frame(type = "ref", year=an[y], value = ref_season[y]))

    }
    #aes(fill=type, colour=type)
    ggplot(df_raw, aes(x=year, y=value,  fill=type)) + geom_bar(stat = "identity", position = "dodge")+
      labs(title = paste0(indicateur, " seasonal aggregation"), x="Year", y= "Value")+
      theme_bw()+
      scale_x_continuous(breaks = scales::pretty_breaks(n = 12))


