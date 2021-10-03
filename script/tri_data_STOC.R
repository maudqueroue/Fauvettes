######                                       ######
###### SCRIPT PREMIER TRI DONNEES PASSEREAUX ######
######                                       ######

# Jeu de donnees entrant : donnees brut MNHN

# Jeu de donnees sortant : histoire de vie des sites
#                          Premier tri donnees MNHN


rm(list=ls())

# Packages necessaires
devtools::install_deps(upgrade="never")

# Load fonctions importantes
devtools::load_all()

# DATA MNHN

data_brut <- read.csv2(here::here('Data',"Data_MNHN.csv")) %>%
  dplyr::select('ACTION','BAGUE','DATE','HEURE','ESPECE','SEXE','AGE',"THEME.SESSION","ID_PROG","LAT","LON","DEPT","LOCALITE")

data <- data_brut

#-- 534407

####### ATTENTION MODIF DATA #####
# -> station 43 et 623 ? lier
# -> station 42 et 285 ? lier
# -> station 200 et 615 ? lier
data$ID_PROG[which(data$ID_PROG==43)]<-623
data$ID_PROG[which(data$ID_PROG==42)]<-285
data$ID_PROG[which(data$ID_PROG==200)]<-615


#1 ON NE GARDE QUE LES DONNES STOC et STOC ROZO
# -----------------------------------------------------------------------

unique(data$THEME.SESSION)

data <- data %>%
  dplyr::filter(data$THEME.SESSION == "STOC" | data$THEME.SESSION == "STOC ROZO")

unique(data$THEME.SESSION)
# -- 503904

# 2 QUELLES ANNEES GARDE T ON ?
#------------------------------------------------------------------------

# On ne garde que les annees apres 2001 pour etre coherent avec les donnees EPS
unique(substr(data$DATE,7,10))

data <- data %>%
  dplyr::filter(substr(data$DATE,7,10) > 2000)

unique(substr(data$DATE,7,10))
#-- 422610

#3. Enlever les dates hors STOC : avril-juillet
#--------------------------------------------------------

unique(substr(data$DATE,4,5))

data <- data %>%
  dplyr::filter(substr(data$DATE,4,5) %in% c("04", "05", "06", "07"))

unique(substr(data$DATE,4,5))
#-- 407645

#4. QUELS SITES CONSERVER ?
#------------------------------------------------------

# enlever les lignes sans nom de site
data <- data %>%
  dplyr::filter(!is.na(ID_PROG))

#-- 407644

# ensemble des identifiants des sites
sites <- sort(unique(data$ID_PROG))

# Le nombre d'annees conservees
years <- seq(min(substr(data$DATE,7,10)),max(substr(data$DATE,7,10)),1)

# Combien de captures par sites par annees
hvie_ID_PROG <- matrix(NA,length(sites),length(years))
rownames(hvie_ID_PROG) <-sites
colnames(hvie_ID_PROG) <-years

for(t in 1:length(years)) {
  for(i in 1:length(sites)) {
    s <- dplyr::filter(data, data$ID_PROG==sites[i]  & substr(data$DATE,7,10)==years[t])
    hvie_ID_PROG[i,t] <- dplyr::n_distinct(s$DATE)
    rm(s)
  }
}

rm(i, t)

hvie_ID_PROG <- dplyr::as_tibble(hvie_ID_PROG) %>%
  tibble::add_column(ID_PROG = sites)

# Sauvegarde des histoire de capture des sites
save(hvie_ID_PROG,file=here::here("output","hvie_ID_PROG.RData"))

# Combien d'annees les sites ont fait au moins 3 captures secondaires ?
hvie_ID_PROG <- hvie_ID_PROG %>%
  tibble::add_column(sum_OC_3 = NA)

for(i in 1:length(sites)) {
  hvie_ID_PROG$sum_OC_3[i] <- length(which(hvie_ID_PROG[i,1:length(years)]>=3))
}

rm(i)

# Quel sont les sites qui n'ont captures (avec 3 occasions secondaires) qu'une annee ou jamais ?

site.rm <- sites[which(hvie_ID_PROG$sum_OC_3<2)]


# On retire les sites qui ne respectent pas le protocole
data <- data %>%
  dplyr::filter(!data$ID_PROG %in% site.rm)
#-- 398523

#5. GARDE QUE LES ESPECES CIBLES
#------------------------------------------------------

# fauvettes
unique(data$ESPECE)

data <- data %>%
  dplyr::filter(data$ESPECE=="SYLATR"|data$ESPECE=="SYLBOR")

unique(data$ESPECE)
#-- 112016

data_STOC <- data

#save file
save(data_STOC,file=here::here("output","data_STOC.RData"))

