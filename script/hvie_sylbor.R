
######                                       ######
###### SCRIPT DONNEES MODELE sylbor          ######
######                                       ######

# Jeu de donnees entrant : histoire de vie des sites
#                          Premier tri donnees MNHN

# Jeu de donnees sortant : hvie_sylbor_tot
#                          hvie_ID_prog_sylbor_tot   


rm(list=ls())

# Packages necessaires
devtools::install_deps(upgrade="never")

# Load fonctions importantes
devtools::load_all() 

# Data
load(here::here("output","data_STOC.RData"))

# 1. Selection de l'espece
#------------------------
# On ne garde que les sylbor
sylbor <- subset(data_STOC,data_STOC$ESPECE=="SYLBOR")
rm(data_STOC)

# 2. Selection de la latitude
#------------------------
CLC_STOC <- read.table(here::here("data","coord_STOC.csv"),head=T,sep=";") %>%
  dplyr::rename(
    long = 'Lon',
    lat  = 'Lat')

ID_PROG <- CLC_STOC %>%
  dplyr::filter(CLC_STOC$lat>45)

ID_to_keep <- unique(ID_PROG$ID_PROG)

sylbor <- sylbor %>%
  dplyr::filter(sylbor$ID_PROG %in% ID_to_keep)

rm(ID_to_keep,CLC_STOC,ID_PROG)

# Informations 
# Annees
years <- seq(min(unique(substr(sylbor$DATE,7,10))),max(unique(substr(sylbor$DATE,7,10))),1)
K <- length(years)
#Identifiant des sylbor 
ID_sylbor <- unique(sylbor$BAGUE)
N <- length(ID_sylbor)

# 3. Redefinition de l'age 
#-----------------------
unique(sylbor$AGE)

sylbor <- Mesanges::new_age(sylbor)

unique(sylbor$new_AGE)

# 4. Creation des histoires de vie
# -----------------------

# hvie
hvie_sylbor <- as.data.frame(matrix(NA,length(ID_sylbor),K))
hvie_sylbor$ID <- ID_sylbor
colnames(hvie_sylbor) <- c(years,"ID")

# Calcul hvie
hvie_sylbor <- Mesanges::hvie_tt(sylbor, hvie_sylbor)
N <- dim(hvie_sylbor)[1]

# On transforme en numerique : Adulte = 2, Jeune = 1
hvie_sylbor[hvie_sylbor=="A"] <- 2
hvie_sylbor[hvie_sylbor=="P"] <- 1
hvie_sylbor[hvie_sylbor=="C"] <- 3
hvie_sylbor[hvie_sylbor=="AP"] <- 4


# 5. On regarde s'il n'y a pas des individus vu jeunes plus tard que la premiere occasion de capture
#----------------------------------

hvie_sylbor <- Mesanges::check_age(hvie_sylbor)
N <- dim(hvie_sylbor)[1]


# 6. On regarde s'il reste des incertains
#------------------------------------------

# si vus une seule fois incertains, on retire la donnees.
check <- which(hvie_sylbor[,1:K]==3) # 18

# Si il est capture qu'une seule annee, on supprime
hvie_sylbor <- Mesanges::check_3(hvie_sylbor, check)

# Est-ce qu'on a fini ? 
check <- which(hvie_sylbor[,1:K]==3) # Oui

rm(check)

# 7. On regarde ceux pour qui on a pas pu trancher entre Adulte et poussin
#----------------------------
check <- which(hvie_sylbor[,1:K]==4) # 43

# On retranscrit en numerique
sylbor$new_AGE[sylbor$new_AGE=="P"]<-1
sylbor$new_AGE[sylbor$new_AGE=="A"]<-2

# S'il ont ete vus qu'une seule annee, on prend le statut defini lors du bagage
# Si pas de baguage la seule annee ou il y a ete vu, on le regarde au cas par cas
# Si individus vus plusieurs annees on regarde au cas par cas
hvie_sylbor <- Mesanges::check_4(sylbor, hvie_sylbor, check)

# Qui doit on trier au cas par cas :
check <- which(hvie_sylbor[,1:K]==4) # 1 individu

#1
ligne <- check[1]%%N
colonne <- ceiling(check[1]/N)
bague <- hvie_sylbor$ID[ligne] #....6562142
subset(sylbor,sylbor$BAGUE==bague)
hvie_sylbor[ligne,]
hvie_sylbor[ligne,colonne] <- 2
hvie_sylbor[ligne,]

# On a fini ?
check <- which(hvie_sylbor[,1:K]==4) #oui

# On supprime les oiseaux qui ne sont finalement pas conserves
# avec uniquement des 0 dans les hvie
hvie_sylbor <- Mesanges::supp_ind(hvie_sylbor)

# On remet à jour N et ID_sylbor
N <- dim(hvie_sylbor)[1]
ID_sylbor <- hvie_sylbor$ID

rm(check,ligne,colonne,bague)

# 8. Gestion des transients
# #-----------------------

hvie_sylbor <- Mesanges::transient_tt(sylbor, hvie_sylbor)
# warnings : pas de probleme

# On supprime les oiseaux qui ne sont finalement pas conserves
# avec uniquement des 0 dans les hvie
hvie_sylbor <- Mesanges::supp_ind(hvie_sylbor)

# On remet à jour N et ID_sylbor
N <- dim(hvie_sylbor)[1]
ID_sylbor <- hvie_sylbor$ID

# 9. Covariable individuelle 
#---------

# Il nous faut le nombre de captures total (-1 pour transience)
# et le nombre d'annees de capture
# pour chaque individus des hvie
hvie_sylbor <- Mesanges::cov_ind_tt(sylbor, hvie_sylbor)


# Quels individus sont problematiques
which(hvie_sylbor$nb_years_capt==0) #0
which(hvie_sylbor$nb_capt==0) #0
which(is.na(hvie_sylbor$nb_years_capt)==T) #0
which(is.na(hvie_sylbor$nb_capt)==T) #0

#Calcul de la covariable :
hvie_sylbor <- Mesanges::calcul_cov_ind(hvie_sylbor)

which(hvie_sylbor$cov_ind<0) # ca parait bon

# 10. Save
#--------------
save(hvie_sylbor,file=here::here('output',"hvie_sylbor.RData"))


# Analyses supplementaires 
#---------------------------
rm(list=ls())
load(here::here('output',"hvie_sylbor_tot_tt.RData"))
load(here::here('output',"hvie_ID_PROG_sylbor_tot_tt.RData"))

info_PROG <- hvie_ID_PROG_sylbor[,-c(1:19)] %>%
  dplyr::mutate(ID_PROG=as.factor(ID_PROG)) 

hvie_sylbor <- hvie_sylbor %>%
  dplyr::mutate(ID_PROG=as.factor(ID_PROG))

rm(hvie_ID_PROG_sylbor)

data <- dplyr::left_join(hvie_sylbor, info_PROG, by = c("ID_PROG")) 

N <- nrow(data)
K <- 19

nb_capt <- NULL
juv <- NULL
ad <- NULL
ad_juv <- NULL
ad_1 <- NULL
juv_1 <- NULL

for(i in 1:N){
  nb_capt[i] <- length(which(data[i,1:K]>0))
  
  juv[i] <- length(which(data[i,1:K]==1))
  ad[i] <- length(which(data[i,1:K]==2))
  
  if (nb_capt[i] >1) {
    if(ad[i]>0 & juv[i]>0){
      ad_juv[i] <- 1 }
    else {ad_juv[i] <- 0}
  }
  
  if(nb_capt[i] == 1) {
    ad_1[i] <- ad[i]
    juv_1[i] <- juv[i]
  }
  else {ad_1[i] <- 0
  juv_1[i] <- 0}
}

length(which(nb_capt>1))
length(which(juv>0))
length(which(ad>0))
length(which(ad_1>0))
length(which(juv_1>0))
length(which(ad_juv==1))

