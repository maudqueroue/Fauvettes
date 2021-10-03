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



