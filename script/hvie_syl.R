######                                       ######
###### SCRIPT DONNEES MODELE PARCAE          ######
######                                       ######

# Jeu de donnees entrant : hvie_sylbor
#                          hvie_sylatr
#                          hvie_ID_PROG

# Jeu de donnees sortant : hvie_sylbor
#                          hvie_sylatr
#                          hvie_ID_syl

rm(list=ls())

# Packages necessaires
devtools::install_deps(upgrade="never")

# Load fonctions importantes
devtools::load_all()

# hvie syl
#--------------------
load(here::here('output',"hvie_sylbor.RData"))
load(here::here('output',"hvie_sylatr.RData"))

# Hvie site
#--------------------
load(here::here('output',"hvie_ID_PROG.RData"))

#Quels sont les sites concernant les parcae
hvie_ID_PROG_sylbor<- hvie_ID_PROG %>%
  dplyr::filter(hvie_ID_PROG$ID_PROG %in% unique(hvie_sylbor$ID_PROG))

#Quels sont les sites concernant les parmaj
hvie_ID_PROG_sylatr <- hvie_ID_PROG %>%
  dplyr::filter(hvie_ID_PROG$ID_PROG %in% unique(hvie_sylatr$ID_PROG))

# Sites en commun pour les deux
#--------------------
hvie_ID_PROG_syl <- hvie_ID_PROG_sylatr %>%
  dplyr::filter(hvie_ID_PROG_sylatr$ID_PROG %in% hvie_ID_PROG_sylbor$ID_PROG)
hvie_ID_PROG_syl <- hvie_ID_PROG_syl[-which(hvie_ID_PROG_syl$ID_PROG==56),]


rm(hvie_ID_PROG, hvie_ID_PROG_sylbor, hvie_ID_PROG_sylatr)

# On filtre pour ne garder qu eles hvie communes aux deux espèces
#--------------------
ID_PROG_syl <- unique(hvie_ID_PROG_syl$ID_PROG)

hvie_sylbor <- hvie_sylbor %>%
  dplyr::filter(hvie_sylbor$ID_PROG %in% ID_PROG_syl)
hvie_sylatr <- hvie_sylatr %>%
  dplyr::filter(hvie_sylatr$ID_PROG %in% ID_PROG_syl)


# 11. Creation vecteur pour lier hvie individus et hvie site
#---------
hvie_sylbor <- Fauvettes::link_hvie(hvie_sylbor, hvie_ID_PROG_syl)
hvie_sylatr <- Fauvettes::link_hvie(hvie_sylatr, hvie_ID_PROG_syl)

save(hvie_ID_PROG_syl,file=here::here('output',"hvie_ID_PROG_syl.RData"))
save(hvie_sylbor,file=here::here('output',"hvie_sylbor.RData"))
save(hvie_sylatr,file=here::here('output',"hvie_sylatr.RData"))


# Infos

data <- hvie_sylbor

N <- unique(data$ID)
K <- 19

nb_capt <- NULL
juv <- NULL
ad <- NULL
ad_rec <- NULL

for(i in 1:nrow(data)){
  nb_capt[i] <- length(which(data[i,1:K]>0))

  juv[i] <- length(which(data[i,1:K]==1))
  ad[i] <- length(which(data[i,1:K]==2))

  if (nb_capt[i] >1) {
    if(ad[i]>0 & juv[i]>0){
      ad_juv[i] <- 1 }
    else {ad_juv[i] <- 0}

    if(ad[i]>1){
      ad_rec[i] <- ad[i] - 1 }
    else {ad_rec[i] <- 0}
  }

}
length(N)
length(which(ad_juv==1))
sum(ad_rec,na.rm=T)

# Infos

data <- hvie_sylatr

N <- unique(data$ID)
K <- 19

nb_capt <- NULL
juv <- NULL
ad <- NULL
ad_rec <- NULL

for(i in 1:nrow(data)){
  nb_capt[i] <- length(which(data[i,1:K]>0))

  juv[i] <- length(which(data[i,1:K]==1))
  ad[i] <- length(which(data[i,1:K]==2))

  if (nb_capt[i] >1) {
    if(ad[i]>0 & juv[i]>0){
      ad_juv[i] <- 1 }
    else {ad_juv[i] <- 0}

    if(ad[i]>1){
      ad_rec[i] <- ad[i] - 1 }
    else {ad_rec[i] <- 0}
  }

}
length(N)
length(which(ad_juv==1))
sum(ad_rec,na.rm=T)


