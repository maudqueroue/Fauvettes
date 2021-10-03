
######                                       ######
###### SCRIPT DONNEES MODELE SYLATR          ######
######                                       ######

# Jeu de donnees entrant : histoire de vie des sites
#                          Premier tri donnees MNHN

# Jeu de donnees sortant : hvie_sylatr_tot
#                          hvie_ID_prog_sylatr_tot   


rm(list=ls())

# Packages necessaires
devtools::install_deps(upgrade="never")

# Load fonctions importantes
devtools::load_all() 

# Data
load(here::here("output","data_STOC.RData"))


# 1. Selection de l'espece
#------------------------
# On ne garde que les sylatr
sylatr <- subset(data_STOC,data_STOC$ESPECE=="SYLATR")
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

sylatr <- sylatr %>%
  dplyr::filter(sylatr$ID_PROG %in% ID_to_keep)

rm(ID_to_keep,CLC_STOC,ID_PROG)

# Informations 
# Annees
years <- seq(min(unique(substr(sylatr$DATE,7,10))),max(unique(substr(sylatr$DATE,7,10))),1)
K <- length(years)
#Identifiant des sylatr 
ID_sylatr <- unique(sylatr$BAGUE)
N <- length(ID_sylatr)

# 3. Redefinition de l'age 
#-----------------------
unique(sylatr$AGE)

sylatr <- Mesanges::new_age(sylatr)

unique(sylatr$new_AGE)

# 4. Creation des histoires de vie
# -----------------------

# hvie
hvie_sylatr <- as.data.frame(matrix(NA,length(ID_sylatr),K))
hvie_sylatr$ID <- ID_sylatr
colnames(hvie_sylatr) <- c(years,"ID")

# Calcul hvie
hvie_sylatr <- Mesanges::hvie_tt(sylatr, hvie_sylatr)
N <- dim(hvie_sylatr)[1]

# On transforme en numerique : Adulte = 2, Jeune = 1
hvie_sylatr[hvie_sylatr=="A"] <- 2
hvie_sylatr[hvie_sylatr=="P"] <- 1
hvie_sylatr[hvie_sylatr=="C"] <- 3
hvie_sylatr[hvie_sylatr=="AP"] <- 4


# 5. On regarde s'il n'y a pas des individus vus jeunes plus tard que la premiere occasion de capture
#----------------------------------

hvie_sylatr <- Mesanges::check_age(hvie_sylatr)
N <- dim(hvie_sylatr)[1]

# 6. On regarde s'il reste des incertains
#------------------------------------------

# si vus une seule fois incertains, on retire la donnees.
check <- which(hvie_sylatr[,1:K]==3) # 38

# Si il est capture qu'une seule annee, on supprime
hvie_sylatr <- Mesanges::check_3(hvie_sylatr, check)

# Qui doit on trier au cas par cas :
check <- which(hvie_sylatr[,1:K]==3) # 3

# cas pas cas
ligne <- check[1]%%N
colonne <- ceiling(check[1]/N)
bague <- hvie_sylatr$ID[ligne] # ....4440015
subset(sylatr,sylatr$BAGUE==bague)
hvie_sylatr[ligne,colonne] <- 1

ligne <- check[2]%%N
colonne <- ceiling(check[2]/N)
bague <- hvie_sylatr$ID[ligne] #....5225935
subset(sylatr,sylatr$BAGUE==bague)
hvie_sylatr[ligne,colonne] <- 1

ligne <- check[3]%%N
colonne <- ceiling(check[3]/N)
bague <- hvie_sylatr$ID[ligne] #....8255330
subset(sylatr,sylatr$BAGUE==bague)
hvie_sylatr[ligne,colonne] <- 1

# Est-ce qu'on a fini ? 
check <- which(hvie_sylatr[,1:K]==3) # Oui

rm(colonne,ligne,bague,check)

# 7. On regarde ceux pour qui on a pas pu trancher entre Adulte et poussin
#----------------------------
check <- which(hvie_sylatr[,1:K]==4) # 110

# On retranscrit en numerique
sylatr$new_AGE[sylatr$new_AGE=="P"]<-1
sylatr$new_AGE[sylatr$new_AGE=="A"]<-2

# S'il ont ete vus qu'une seule annee, on prend le statut defini lors du bagage
# Si pas de baguage la seule annee ou il y a ete vu, on le regarde au cas par cas
# Si individus vus plusieurs annees on regarde au cas par cas
hvie_sylatr <- Mesanges::check_4(sylatr, hvie_sylatr, check)

# Qui doit on trier au cas par cas :
check <- which(hvie_sylatr[,1:K]==4) # 2 individu

#1
ligne <- check[1]%%N
colonne <- ceiling(check[1]/N)
bague <- hvie_sylatr$ID[ligne] #....5682208
subset(sylatr,sylatr$BAGUE==bague)
hvie_sylatr[ligne,]
hvie_sylatr[ligne,colonne] <- 2
hvie_sylatr[ligne,]

#2
ligne <- check[2]%%N
colonne <- ceiling(check[2]/N)
bague <- hvie_sylatr$ID[ligne] #....5566918
subset(sylatr,sylatr$BAGUE==bague)
hvie_sylatr[ligne,]
hvie_sylatr[ligne,colonne] <- 1
hvie_sylatr[ligne,]

# On a fini ?
check <- which(hvie_sylatr[,1:K]==4) #oui

# On supprime les oiseaux qui ne sont finalement pas conserves
# avec uniquement des 0 dans les hvie
hvie_sylatr <- Mesanges::supp_ind(hvie_sylatr)

# On remet à jour N et ID_sylatr
N <- dim(hvie_sylatr)[1]
ID_sylatr <- hvie_sylatr$ID

rm(check,ligne,colonne,bague)

# 8. Gestion des transients
# #-----------------------

hvie_sylatr <- Mesanges::transient_tt(sylatr, hvie_sylatr)
# warnings : pas de probleme

# On supprime les oiseaux qui ne sont finalement pas conserves
# avec uniquement des 0 dans les hvie
hvie_sylatr <- Mesanges::supp_ind(hvie_sylatr)

# On remet à jour N et ID_sylatr
N <- dim(hvie_sylatr)[1]
ID_sylatr <- hvie_sylatr$ID

# 9. Covariable individuelle 
#---------

# Il nous faut le nombre de captures total (-1 pour transience)
# et le nombre d'annees de capture
# pour chaque individus des hvie
hvie_sylatr <- Mesanges::cov_ind_tt(sylatr, hvie_sylatr)


# Quels individus sont problematiques
which(hvie_sylatr$nb_years_capt==0) #0
which(hvie_sylatr$nb_capt==0) #0
which(is.na(hvie_sylatr$nb_years_capt)==T) #0
which(is.na(hvie_sylatr$nb_capt)==T) #0

#Calcul de la covariable :
hvie_sylatr <- Mesanges::calcul_cov_ind(hvie_sylatr)

which(hvie_sylatr$cov_ind<0) # ca parait bon

# 10. Save
#--------------
save(hvie_sylatr,file=here::here('output',"hvie_sylatr.RData"))


# Analyses supplementaires 
#---------------------------
rm(list=ls())
load(here::here('output',"hvie_sylatr_tot_tt.RData"))
load(here::here('output',"hvie_ID_PROG_sylatr_tot_tt.RData"))

info_PROG <- hvie_ID_PROG_sylatr[,-c(1:19)] %>%
  dplyr::mutate(ID_PROG=as.factor(ID_PROG)) 

hvie_sylatr <- hvie_sylatr %>%
  dplyr::mutate(ID_PROG=as.factor(ID_PROG))

rm(hvie_ID_PROG_sylatr)

data <- dplyr::left_join(hvie_sylatr, info_PROG, by = c("ID_PROG")) 

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

