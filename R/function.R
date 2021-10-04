#' give point
#'
#' @param  ID
#'
#' @return
#' @export
#'
give_point <- function(ID) {
  dsf <- ID %>%
    dplyr::select("long","lat") %>%
    sf::st_as_sf(coords = c("long","lat"), crs = 4326) %>%
    sf::st_transform(crs = 2154)
}

#' Give point by site
#'
#' @param  dsf_STOC, dsf_EPS, dist_km
#'
#' @return
#' @export
#'
give_point_by_site <- function(dsf_STOC_i, dist_km, CLC_EPS, dsf_EPS) {
  coord_buffer_STOC <- sf::st_buffer(dsf_STOC_i, dist=units::set_units(dist_km, km))
  int <- sf::st_intersects(coord_buffer_STOC$geometry, dsf_EPS)
  pts <- CLC_EPS[int[[1]],"point"]
  return(pts)
}

#' @param  data
#'
#' @return
#' @export
#'
supr_point <- function (data, sp) {

  if(sp =="SYL") {
    # Quels points n ont jamais contacte l'sp ?
    sum_sp_point <- data %>%
      dplyr::group_by(point) %>%
      dplyr::summarise_at(.vars= dplyr::vars("SYLATR","SYLBOR"), .funs= "sum")
    point_supr <- sum_sp_point$point[which(sum_sp_point$SYLATR == 0 | sum_sp_point$SYLBOR ==0)]
  }

  if(sp =="SYLBOR") {
    # Quels points n ont jamais contacte l'sp ?
    sum_sp_point <- data %>%
      dplyr::group_by(point) %>%
      dplyr::summarise(n= sum(SYLBOR))
    point_supr <- sum_sp_point$point[which(sum_sp_point$n == 0)]
  }

  if(sp =="SYLATR") {
    # Quels points n ont jamais contacte l'sp ?
    sum_sp_point <- data %>%
      dplyr::group_by(point) %>%
      dplyr::summarise(n= sum(SYLATR))
    point_supr <- sum_sp_point$point[which(sum_sp_point$n == 0)]
  }

  # On supprime ces points
  data <- data %>%
    dplyr::filter(!(data$point %in% point_supr))

  return(data)

}


#' Nouveau calcul d index
#'
#' @param  data
#'
#' @return
#' @export
#'
index_diff <- function(data, sp) {


  years <- seq(min(data$annee),max(data$annee),1)
  K <- length(years)
  ID_point <- unique(data$point)
  N <- length(ID_point)
  hvie_point <- as.data.frame(matrix(NA,N,K))
  hvie_point$ID <- ID_point


  for (i in 1:N) {
    data_point <- subset(data,data$point==hvie_point$ID[i])
    for (j in 1:K) {
      data_point_annee <- subset(data_point,data_point$annee==years[j])

      if(sp == "SYLATR") {

        if (nrow(data_point_annee)==0) {hvie_point[i,j] <- NA}
        if (nrow(data_point_annee)>=1) {hvie_point[i,j] <- data_point_annee$SYLATR[1]}
        rm(data_point_annee)
      }

      if(sp == "SYLBOR") {

        if (nrow(data_point_annee)==0) {hvie_point[i,j] <- NA}
        if (nrow(data_point_annee)>=1) {hvie_point[i,j] <- data_point_annee$SYLBOR[1]}
        rm(data_point_annee)
      }
    }

    rm(data_point)
  }


  fluct <- NULL

  for (j in 2:K) {
    sub_data <- subset(hvie_point,is.na(hvie_point[,(j-1)])==F & is.na(hvie_point[,j])==F)
    fluct[j] <- sum(sub_data[,j]) / sum(sub_data[,(j-1)])
  }

  ens <- NULL
  ens[1] <- 1
  for (i in 2:K) {
    ens[i] <- ens[i-1] * fluct[i]
  }

  return(ens)
}


#' Nb individus
#'
#' @param  data
#'
#' @return
#' @export
#'
calcul_n <- function(data,sp,index) {

  years <- seq(min(data$annee),max(data$annee),1)
  K <- length(years)
  ID_point <- unique(data$point)
  N <- length(ID_point)
  hvie_point <- as.data.frame(matrix(NA,N,K))
  hvie_point$ID <- ID_point

  for (i in 1:N) {
    data_point <- subset(data,data$point==hvie_point$ID[i])
    for (j in 1:K) {
      data_point_annee <- subset(data_point,data_point$annee==years[j])

      if(sp == "SYLATR") {

        if (nrow(data_point_annee)==0) {hvie_point[i,j] <- NA}
        if (nrow(data_point_annee)>=1) {hvie_point[i,j] <- data_point_annee$SYLATR[1]}
        rm(data_point_annee)
      }

      if(sp == "SYLBOR") {

        if (nrow(data_point_annee)==0) {hvie_point[i,j] <- NA}
        if (nrow(data_point_annee)>=1) {hvie_point[i,j] <- data_point_annee$SYLBOR[1]}
        rm(data_point_annee)
      }
    }
    rm(data_point)
  }

  nb_ind_an <- apply(hvie_point[,1:19],2,sum,na.rm=T)

  nb_point_an <- NULL
  for(i in 1:K) {
    nb_point_an[i] <- length(which(is.na(hvie_point[,i])==F))
  }

  n_max <- which(nb_point_an==max(nb_point_an))
  nm <- nb_ind_an[n_max]/index[n_max]

  return(round(as.numeric(nm)))
}


#' Calcul index
#'
#' @param  data
#'
#' @return
#' @export
#'
index_new <- function(data, sp) {

  # Nombre de stations STOC
  nsites <- length(unique(data$point))
  # Nombre d'annees
  nyears <- length(unique(data$annee))

  # effet fixe annee et site

  if(sp == "SYLATR"){
    model.glm <- glm(as.numeric(SYLATR) ~ 0 + as.factor(point) + as.factor(annee), data = data, family = 'poisson')
  }

  if(sp == "SYLBOR"){
    model.glm <- glm(as.numeric(SYLBOR) ~ 0 + as.factor(point) + as.factor(annee), data = data, family = 'poisson')
  }

  # Calcul de l'index en supposant que alpha et beta suivent des normales et en faisant du Monte Carlo, on simule un grand nombre de valeurs.

  alphaihat <- model.glm$coefficients[1:nsites]
  se_alphaihat <- summary(model.glm)$coefficients[,'Std. Error'][1:nsites]

  nbMC <- 100
  aik <- matrix(NA, nrow = nsites, ncol = nbMC)
  for (i in 1:nsites){
    for (j in 1:nbMC){
      aik[i,j] <- rnorm(1, mean = alphaihat[i], sd = se_alphaihat[i])
    }
  }

  betathat <- model.glm$coefficients[(nsites+1):length(model.glm$coefficients)]
  se_betathat <- summary(model.glm)$coefficients[,'Std. Error'][(nsites+1):length(model.glm$coefficients)]

  btk <- matrix(NA, nrow = nyears - 1, ncol = nbMC)
  for (i in 1:(nyears - 1)){
    for (j in 1:nbMC){
      btk[i,j] <- rnorm(1, mean = betathat[i], sd = se_betathat[i])
    }
  }

  #gammaitk
  gammaitk <- array(NA, dim = c(nsites, nyears - 1, nbMC))
  for (t in 1:(nyears - 1)){
    for (i in 1:nsites){
      for (k in 1:nbMC){
        gammaitk[i,t,k] <- exp(btk[t,k] + aik[i,k])
      }

    }
  }

  #gammat
  gammat <- t(apply(gammaitk,c(2,3), sum, na.rm = TRUE))

  # On calcule la moyenne des log et leur variance
  logyt <- apply(log(gammat), 2, mean)

  sigma2tt <- (log(gammat) - matrix(rep(logyt, nbMC), nrow = nbMC, byrow = T))^2
  sigma2t <- apply(sigma2tt, 2, mean)


  index_mean <- logyt
  index_var <- sigma2t

  out <- list(index_mean, index_var, gammaitk)
  return(out)
}


#' Calcul index
#'
#' @param  data
#'
#' @return
#' @export
#'
index_new_carre <- function(data_n, sp) {


  # effet fixe annee et site

  if(sp == "SYLATR"){
    data <- data_n %>%
      dplyr::group_by(carre,annee) %>%
      dplyr::summarise(n_ind = sum(SYLATR))
  }

  if(sp == "SYLBOR"){
    data <- data_n %>%
      dplyr::group_by(carre,annee) %>%
      dplyr::summarise(n_ind = sum(SYLBOR))

  }

  # Nombre de stations STOC
  nsites <- length(unique(data$carre))
  # Nombre d'annees
  nyears <- length(unique(data$annee))

  model.glm <- glm(as.numeric(n_ind) ~ 0 + as.factor(carre) + as.factor(annee), data = data, family = 'poisson')


  # Calcul de l'index en supposant que alpha et beta suivent des normales et en faisant du Monte Carlo, on simule un grand nombre de valeurs.

  alphaihat <- model.glm$coefficients[1:nsites]
  se_alphaihat <- summary(model.glm)$coefficients[,'Std. Error'][1:nsites]

  nbMC <- 100
  aik <- matrix(NA, nrow = nsites, ncol = nbMC)
  for (i in 1:nsites){
    for (j in 1:nbMC){
      aik[i,j] <- rnorm(1, mean = alphaihat[i], sd = se_alphaihat[i])
    }
  }

  betathat <- model.glm$coefficients[(nsites+1):length(model.glm$coefficients)]
  se_betathat <- summary(model.glm)$coefficients[,'Std. Error'][(nsites+1):length(model.glm$coefficients)]

  btk <- matrix(NA, nrow = nyears - 1, ncol = nbMC)
  for (i in 1:(nyears - 1)){
    for (j in 1:nbMC){
      btk[i,j] <- rnorm(1, mean = betathat[i], sd = se_betathat[i])
    }
  }

  #gammaitk
  gammaitk <- array(NA, dim = c(nsites, nyears - 1, nbMC))
  for (t in 1:(nyears - 1)){
    for (i in 1:nsites){
      for (k in 1:nbMC){
        gammaitk[i,t,k] <- exp(btk[t,k] + aik[i,k])
      }

    }
  }

  #gammat
  gammat <- t(apply(gammaitk,c(2,3), sum, na.rm = TRUE))

  # On calcule la moyenne des log et leur variance
  logyt <- apply(log(gammat), 2, mean)

  sigma2tt <- (log(gammat) - matrix(rep(logyt, nbMC), nrow = nbMC, byrow = T))^2
  sigma2t <- apply(sigma2tt, 2, mean)


  index_mean <- logyt
  index_var <- sigma2t

  out <- list(index_mean, index_var)
  return(out)
}

#' Calcul nouvel age
#'
#' @param  data
#'
#' @return
#' @export
#'
new_age <- function (data) {

  data <- data %>%
    tibble::add_column(new_AGE = NA)

  # Juv
  data$new_AGE[which(data$AGE=="PUL")]<- "P"
  data$new_AGE[which(data$AGE=="1A")] <- "P"
  data$new_AGE[which(data$AGE=="1A?")]<- "P"

  # Ad
  data$new_AGE[which(data$AGE=="+1A")]<- "A"
  data$new_AGE[which(data$AGE=="+1?")]<- "A"
  data$new_AGE[which(data$AGE=="2A")] <- "A"
  data$new_AGE[which(data$AGE=="2A?")]<- "A"
  data$new_AGE[which(data$AGE=="+2A")]<- "A"
  data$new_AGE[which(data$AGE=="+2?")]<- "A"
  data$new_AGE[which(data$AGE=="3A")] <- "A"
  data$new_AGE[which(data$AGE=="+3A")]<- "A"
  data$new_AGE[which(data$AGE=="4A")] <- "A"
  data$new_AGE[which(data$AGE=="+4A")]<- "A"
  data$new_AGE[which(data$AGE=="+6A")]<- "A"

  # Incertains
  data$new_AGE[which(data$AGE=="")]   <- "C"
  data$new_AGE[which(data$AGE=="VOL")]<- "C"

  return(data)
}

#' histoires de capture
#'
#' @param  data
#'
#' @return
#' @export
#'
hvie <- function(data, hvie_sp) {

  K <- 19
  N <- nrow(hvie_sp)
  hvie_sp <- hvie_sp %>%
    tibble::add_column(censure = 1) %>%
    tibble::add_column(ID_PROG = NA)


  for (i in 1:N) {
    data_ind <- subset(data,data$BAGUE==hvie_sp$ID[i])
    nb_site <- unique(data_ind$ID_PROG)

    a <- subset(data_ind,data_ind$ID_PROG==nb_site[1])
    hvie_sp$ID_PROG[i] <- nb_site[1]

    for (j in 1:K) {
      b <- subset(a,substr(a$DATE,7,10)==years[j])

      if (nrow(b)==0) {hvie_sp[i,j] <- 0}
      if (nrow(b)==1) {hvie_sp[i,j] <- b$new_AGE[1]}
      if (nrow(b)>1)  {

        c <- unique(b$new_AGE)

        if(length(which(c=="A"))>0  & length(which(c=="P"))>0 ) {hvie_sp[i,j]  <- "AP"}
        if(length(which(c=="A"))>0  & length(which(c=="P"))==0) {hvie_sp[i,j]  <- "A"  }
        if(length(which(c=="A"))==0 & length(which(c=="P"))>0 ) {hvie_sp[i,j]  <- "P"  }
        if(length(which(c=="A"))==0 & length(which(c=="P"))==0 ){hvie_sp[i,j]  <- "C"  }

        rm(c)
      }
      rm(b)
    }
    rm(a)

    if (length(nb_site)>1) {

      for (s in 2:length(nb_site)) {

        hvie_new <- c(rep(0,K),data_ind$BAGUE[1],1,nb_site[s])
        a <- subset(data_ind,data_ind$ID_PROG==nb_site[s])

        for (j in 1:K) {
          b <- subset(a,substr(a$DATE,7,10)==years[j])

          if (nrow(b)==0) {hvie_new[j] <- 0}
          if (nrow(b)==1) {hvie_new[j] <- b$new_AGE[1]}
          if (nrow(b)>1)  {

            c <- unique(b$new_AGE)

            if(length(which(c=="A"))>0  & length(which(c=="P"))>0 ) {hvie_new[j]  <- "AP"}
            if(length(which(c=="A"))>0  & length(which(c=="P"))==0) {hvie_new[j]  <- "A" }
            if(length(which(c=="A"))==0 & length(which(c=="P"))>0 ) {hvie_new[j]  <- "P" }
            if(length(which(c=="A"))==0 & length(which(c=="P"))==0 ){hvie_new[j]  <- "C" }
            rm(c)
          }
          rm(b)
        }
        rm(a)

        hvie_sp <- rbind(hvie_sp, hvie_new)
        rm(hvie_new)
      }
    }
  }

  return(hvie_sp)

}


#' Calcul transient
#'
#' @param  data, hvie_sp
#'
#' @return
#' @export
#'
transient <- function(data, hvie_sp) {

  K <- 19
  N <- nrow(hvie_sp)

  # calcul nombre de captures
  transient <- matrix(NA,N,K)

  for (i in 1:N) {
    a <- subset(data,data$BAGUE==hvie_sp$ID[i] & data$ID_PROG==hvie_sp$ID_PROG[i])

    for (j in 1:K) {
      b <- subset(a,substr(a$DATE,7,10)==years[j])
      transient[i,j] <- nrow(b)
      rm(b)
    }
    rm(a)
  }

  rm(i,j)


  # calcul des nb de captures
  juv <- NULL
  ad <- NULL
  ad_juv <- NULL

  for(i in 1:N){

    juv[i] <- length(which(hvie_sp[i,1:K]==1))
    ad[i] <- length(which(hvie_sp[i,1:K]==2))

    if(ad[i]>0 & juv[i]>0){ad_juv[i] <- 1 }
    else{ad_juv[i] <- 0}
  }


  # Première occasion de capture adulte
  e <- NULL
  e_ad <- NULL
  e_juv <- NULL

  for (i in 1:N){
    temp <- 1:K
    e <- c(e,min(temp[hvie_sp[i,1:K]>=1]))
    e_ad <- c(e_ad,min(temp[hvie_sp[i,1:K]==2]))
    e_juv <- c(e_juv,min(temp[hvie_sp[i,1:K]==1]))}

  # Gestion transients
  for (i in 1:N){
    if(ad_juv[i] == 0) {
      if(transient[i,e[i]]==1) {
        if (hvie_sp[i,e[i]]==1) {next}
        else {hvie_sp[i,e[i]]<-0}
      }

      if(transient[i,e[i]]>1) {next}
    }

    if(ad_juv[i] == 1) {

      hvie_new <- hvie_sp[i,]
      hvie_new$censure <- -1

      hvie_sp[i,e_juv[i]] <- 0

      if(e_ad[i] < 19) {
        hvie_new[,(e_ad[i]+1):19] <- rep(0,(19-e_ad[i]))
      }

      if(e_ad[i] == 19) {hvie_new <- hvie_new}

      hvie_sp <- rbind(hvie_sp, hvie_new)
      rm(hvie_new)

    }
  }

  for(i in which(ad_juv==1)){
    if(transient[i,e_ad[i]]==1) {hvie_sp[i,e_ad[i]]<-0}
    if(transient[i,e_ad[i]]>1) {next}
  }

  return(hvie_sp)
}


#' Suppression individu hvie = 0
#'
#' @param  data
#'
#' @return
#' @export
#'
supp_ind <- function(hvie_sp) {
  # On verifie le nombre d'individus
  N <- dim(hvie_sp)[1]

  # On supprime les oiseaux qui ne sont vu que transients
  hvie_sp$sum_hvie <- rep(NA,N)
  for(i in 1:N){
    hvie_sp$sum_hvie[i] <-sum(as.numeric(hvie_sp[i,1:K]))
  }

  # On supprime les choses inutiles
  hvie_sp <- hvie_sp %>%
    dplyr::filter(!hvie_sp$sum_hvie==0) %>%
    dplyr::select(!sum_hvie)

  return(hvie_sp)
}


#' Gestion de l'age : forcement adulte apres premiere occasion de capture
#'
#' @param  hvie_sp
#'
#' @return
#' @export
#'
check_age <- function(hvie_sp) {
  K <- 19
  e <- NULL
  N <- nrow(hvie_sp)

  for (i in 1:N){
    temp <- 1:K
    e <- c(e,min(temp[hvie_sp[i,1:K]>=1]))
  }

  for (i in 1:N){

    if(e[i] < 19){
      for (j in (e[i]+1):K){
        if (hvie_sp[i,j] == 1) {hvie_sp[i,j] <- 2}
        if (hvie_sp[i,j] == 3) {hvie_sp[i,j] <- 2}
        if (hvie_sp[i,j] == 4) {hvie_sp[i,j] <- 2}
        if (hvie_sp[i,j] == 2) {next}
        if (hvie_sp[i,j] == 0) {next}
      }
    }
    if(e[i]==19) {next}
  }
  return(hvie_sp)
}


#' Gestion des incertains : Si il est capture qu'une seule annee, on supprime
#'
#' @param  data
#'
#' @return
#' @export
#'
check_3 <- function(hvie_sp, check) {
  N <- dim(hvie_sp)[1]
  for(i in 1:length(check)){
    ligne <- check[i]%%N
    colonne <- ceiling(check[i]/N)
    if(length(which(as.numeric(hvie_sp[ligne,1:K])>0))==1){hvie_sp[ligne,colonne] <- 0}
    else{next}
    rm(ligne,colonne)
  }
  return(hvie_sp)
}


#' Gestion quand plusieurs age la même annee
#'
#' @param  data
#'
#' @return
#' @export
#'
check_4 <- function(data, hvie_sp, check) {

  # S'ils ont ete vus qu'une seule annee, on prend le statut defini lors du bagage
  # Si pas de baguage la seule annee, on regarde au cas par cas

  N <- dim(hvie_sp)[1]

  for(i in 1:length(check)){
    ligne <- check[i]%%N
    colonne <- ceiling(check[i]/N)
    bague <- hvie_sp$ID[ligne]
    a <- subset(data,data$BAGUE==bague)
    b <- subset(a, a$ACTION=="B")
    if(nrow(b)>0) {hvie_sp[ligne,colonne] <-b$new_AGE[1]}
    else {next}
    rm(ligne, colonne, bague, a, b)
  }
  return(hvie_sp)
}


#' a quel site appartient chaque oiseau
#'
#' @param  data
#'
#' @return
#' @export
#'
ind_site <- function(data, hvie_sp) {

  N <- dim(hvie_sp)[1]
  hvie_sp$ID_PROG <- rep(NA,N)
  hvie_sp$nb_site <- rep(NA,N)

  for (i in 1:N){
    # On ne garde dans les donnees correspondant a l'oiseau i
    a <- subset(data, data$BAGUE==hvie_sp$ID[i])
    # On r?cup?re les ou l'identifiant(s) de la station
    hvie_sp$ID_PROG[i] <- paste(unique(a$ID_PROG),collapse="/")
    #On regarde s'il y a differents station pour un meme oiseaux
    hvie_sp$nb_site[i] <- length(unique(a$ID_PROG))
    rm(a)
  }

  return(hvie_sp)
}

#' covariable individuelle
#'
#' @param  data
#'
#' @return
#' @export
#'
cov_ind <- function(data, hvie_sp) {

  # Il nous faut le nombre de captures total (-1 pour transience si premiere capt adulte)
  # et le nombre d'annees de capture pour chaque individus des hvie

  N <- dim(hvie_sp)[1]
  K <- 19
  hvie_sp$nb_capt <- rep(NA,N)
  hvie_sp$nb_years_capt <- rep(NA,N)

  # calcul nombre de captures
  nb_capt <- matrix(NA,N,K)

  for (i in 1:N) {
    a <- subset(data,data$BAGUE==hvie_sp$ID[i] & data$ID_PROG==hvie_sp$ID_PROG[i])

    for (j in 1:K) {
      b <- subset(a,substr(a$DATE,7,10)==years[j])
      nb_capt[i,j] <- nrow(b)
      rm(b)
    }
    rm(a)
  }

  rm(i,j)

  e <- NULL
  for (i in 1:N){
    temp <- 1:K
    e <- c(e,min(temp[hvie_sp[i,1:K]>=1]))
  }

  for (i in 1:N) {
    hvie_sp$nb_years_capt[i] <- length(which(as.numeric(hvie_sp[i,1:K])>0))
    if(hvie_sp[i,e[i]]==2 & nb_capt[i,e[i]]>1){
      hvie_sp$nb_capt[i] <- sum(nb_capt[i,which(as.numeric(hvie_sp[i,1:K])>0)])-1}
    else{
      hvie_sp$nb_capt[i] <- sum(nb_capt[i,which(as.numeric(hvie_sp[i,1:K])>0)])}
  }
  return(hvie_sp)
}

#' calcul covariable individuelle
#'
#' @param  data
#'
#' @return
#' @export
#'
calcul_cov_ind <- function(hvie_sp) {

  N <- dim(hvie_sp)[1]
  hvie_sp$cov_ind <- rep(NA,N)

  for(i in 1:N) {
    hvie_sp$cov_ind[i] <- log(1+((hvie_sp$nb_capt[i]-hvie_sp$nb_years_capt[i])/hvie_sp$nb_years_capt[i]))
  }
  hvie_sp <- hvie_sp %>%
    dplyr::select(!nb_capt) %>%
    dplyr::select(!nb_years_capt)
  return(hvie_sp)
}

#' lien hvie_PROG et hvie_IND
#'
#' @param  data
#'
#' @return
#' @export
#'
link_hvie <- function(hvie_sp, hvie_PROG) {

  # On cree le vecteur avec le nouveau numero de site de 1 au nombre total de sites
  N <- dim(hvie_sp)[1]

  hvie_sp$vector_ID_PROG <- rep(NA,N)

  for(i in 1:N){
    hvie_sp$vector_ID_PROG[i] <- which(hvie_PROG$ID_PROG==hvie_sp$ID_PROG[i])
  }
  return(hvie_sp)
}

