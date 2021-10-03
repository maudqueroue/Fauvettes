rm(list=ls())

# Packages necessaires
devtools::install_deps(upgrade="never")

# Load fonctions importantes
devtools::load_all() 


# Carte France
#-------------------

Fr <- sf::st_read(dsn = here::here('data','France'), 
                  layer = "DEPARTEMENT",
                  quiet = TRUE) %>% 
  tibble::add_column(PAYS = "FR") %>% 
  sf::st_transform(4326) %>% 
  dplyr::group_by(PAYS)%>% 
  dplyr::summarize()

# DATA MNHN

data_brut <- read.csv2(here::here('Data',"Data_MNHN.csv")) %>%
  dplyr::select('ACTION','BAGUE','DATE','HEURE','ESPECE','SEXE','AGE',"THEME.SESSION","ID_PROG","LAT","LON","DEPT","LOCALITE")

data <- data_brut 

rm(data_brut)
#-- 534407 

####### ATTENTION MODIF DATA #####
# -> station 43 et 623 ? lier
# -> station 42 et 285 ? lier 
# -> station 200 et 615 ? lier
data$ID_PROG[which(data$ID_PROG==43)]<-623
data$ID_PROG[which(data$ID_PROG==42)]<-285
data$ID_PROG[which(data$ID_PROG==200)]<-615

# on supprime data sans date, sans ID_PROG,sas espece
data <- data %>% 
  dplyr::filter(!is.na(ID_PROG))%>% 
  dplyr::filter(!is.na(DATE)) %>% 
  dplyr::filter(!is.na(ESPECE))


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

data$year <- substr(data$DATE,7,10)

# 2 QUELS MOIS GARDE T ON ?
#------------------------------------------------------------------------
unique(substr(data$DATE,4,5))

data <- data %>%
  dplyr::filter(substr(data$DATE,4,5) %in% c("04", "05", "06", "07")) 

unique(substr(data$DATE,4,5))

# On reformate les données
#------------------------------------------------------------------------
data_STOC <- data %>% 
  dplyr::group_by(ID_PROG,DATE,ESPECE,year) %>%
  dplyr::summarise(n=dplyr::n()) %>%
  tidyr::pivot_wider(names_from = ESPECE, values_from = n) 

# on range dans l'ordre croissant des dates
data_STOC$DATE <- lubridate::dmy(data_STOC$DATE)
data_STOC  <- data_STOC %>%
  dplyr::arrange(DATE)

data_STOC[is.na(data_STOC)] <- 0

# On garde que les sylbor
#------------------------------------------------------------------------
sylbor_STOC <- data_STOC %>%
  dplyr::select('ID_PROG','DATE','year','SYLBOR')


# Combien de sylbor par site par an ? 
#------------------------------------------------------------------------
sylbor_STOC <- sylbor_STOC %>%
  dplyr::group_by(ID_PROG,year) %>%
  dplyr::summarise(n=sum(SYLBOR))

sylbor_STOC <- sylbor_STOC %>%
  dplyr::group_by(ID_PROG) %>%
  dplyr::summarise(n=mean(n))

# Les differents points avec leur longitude et latitude et type CLC
CLC_STOC <- read.table(here::here("output","CLC_STOC.csv"),head=T,sep=";") %>%
  dplyr::select('ID_PROG',"lat",'long')

# On relie les stations avec leur coordonnées
CLC_STOC <- dplyr::left_join(sylbor_STOC, CLC_STOC, by = c("ID_PROG")) %>%
  dplyr::filter(!is.na(long)) %>%
  dplyr::filter(!is.na(lat))

# Creation de la couche de points
dsf_STOC <-  CLC_STOC %>%
  sf::st_as_sf(coords = c("long","lat"), crs = 4326) %>%
  sf::st_transform(crs = 4326)

dsf_STOC_P <- dsf_STOC[which(dsf_STOC$n>0),]

# Plot carte
ggplot2::ggplot(Fr) +
  ggplot2::geom_sf() +
  ggplot2::geom_sf(data = dsf_STOC_P, ggplot2::aes(color=n)) +
  ggplot2::coord_sf(crs = 4326, datum = sf::st_crs(4326)) +
  ggplot2::guides(fill = FALSE) +
  ggplot2::ggtitle("France") +
  ggplot2::theme(title = ggplot2::element_text(size = 16))+
  ggplot2::scale_color_gradientn(colours = rainbow(6))


# Data EPS
#------------------

load(here::here("output","data_EPS.RData"))

data_EPS <- data_EPS[,c('annee','carre','point','SYLBOR')]

# on supprime data sans date, sans point
data_EPS <- data_EPS %>% 
  dplyr::filter(!is.na(point))%>% 
  dplyr::filter(!is.na(annee))

data_EPS[is.na(data_EPS)] <- 0

lien_carre_point <- data_EPS %>%
  dplyr::distinct(carre,.keep_all = TRUE) 

lien_carre_point <- lien_carre_point[,c('carre','point')]

sylbor_EPS <- data_EPS %>%
  dplyr::group_by(annee,carre)  %>%
  dplyr::summarise(n=sum(SYLBOR))

sylbor_EPS <- sylbor_EPS %>%
  dplyr::group_by(carre) %>%
  dplyr::summarise(n=mean(n))

# Les differents points avec leur longitude et latitude et type CLC
load(here::here("output","CLC_EPS.RData"))

# On eneleve les point sans lat ou long
CLC_EPS <- CLC_EPS %>%
  dplyr::select("point", "long", "lat") %>%
  dplyr::filter(!is.na(long)) %>%
  dplyr::filter(!is.na(lat)) 
     
# On relie les point avec leur coordonnées
sylbor_EPS <- dplyr::left_join(sylbor_EPS, lien_carre_point, by = c("carre"))
  
CLC_EPS <- dplyr::left_join(sylbor_EPS, CLC_EPS, by = c("point")) %>%
  dplyr::filter(!is.na(long)) %>%
  dplyr::filter(!is.na(lat)) 

# Creation de la couche de points
dsf_EPS <-  CLC_EPS %>%
  sf::st_as_sf(coords = c("long","lat"), crs = 4326) %>%
  sf::st_transform(crs = 4326)

dsf_EPS_P <- dsf_EPS[which(dsf_EPS$n>0),]

# Plot EPS
ggplot2::ggplot(Fr) +
  ggplot2::geom_sf() +
  ggplot2::geom_sf(data = dsf_EPS_P, ggplot2::aes(color=n)) +
  ggplot2::coord_sf(crs = 4326, datum = sf::st_crs(4326)) +
  ggplot2::guides(fill = FALSE) +
  ggplot2::ggtitle("France") +
  ggplot2::theme(title = ggplot2::element_text(size = 16))+
  ggplot2::scale_color_gradientn(colours = rainbow(6))

# Gam latitude
#----------------------------------------------
library(gamm4)

model <- gamm4(round(n)~s(lat),family=poisson,data=CLC_EPS)

plot.gam(model$gam)
