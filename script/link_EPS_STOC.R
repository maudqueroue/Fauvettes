rm(list=ls())

# Packages necessaires
devtools::install_deps(upgrade="never")

# Load fonctions importantes
devtools::load_all()

# Les differents points EPS avec leur longitude et latitude et type CLC
load(here::here("output","data_EPS.RData"))
geo_EPS <- data_EPS %>%
  dplyr::rename(
    long = 'longitude_wgs84',
    lat  = 'latitude_wgs84') %>%
  dplyr::distinct(point, .keep_all=T) %>%
  dplyr::select("point", "long", "lat") %>%
  dplyr::filter(!is.na(long)) %>%
  dplyr::filter(!is.na(lat))

rm(data_EPS)

# Les differents points STOC avec leur longitude et latitude et type CLC

geo_STOC <- read.table(here::here("data","coord_STOC.csv"),head=T,sep=";") %>%
  dplyr::rename(
    long = 'Lon',
    lat  = 'Lat')


# Creation des couche de points
dsf_STOC <-  Fauvettes::give_point(geo_STOC)
dsf_EPS <-  Fauvettes::give_point(geo_EPS)


# Point par par station STOC
EPS_by_STOC <- list()
for (i in 1:nrow(geo_STOC)) {
  EPS_by_STOC[[i]] <- Fauvettes::give_point_by_site(dsf_STOC[i,], 25, geo_EPS, dsf_EPS)
}

save(EPS_by_STOC, file  = here::here("output","EPS_by_STOC.RData"))


##### Info en supplément

EPS_by_STOC <- list()
for (i in 1:nrow(geo_STOC)) {
  EPS_by_STOC[[i]] <- Fauvettes::give_point_by_site(dsf_STOC[i,], 25, geo_EPS)
}

a <- NULL
for (i in 1:nrow(geo_STOC)) {
  a[i] <- length(EPS_by_STOC[[i]])
}

hist(a, breaks= 50, col= "ivory3", xlab="Nombre de point d'écoute", ylab="Fréquence", main = "Nombre de points d'écoute dans un buffer de 20 km autour de la station STOC")
abline(v=20, lty=2, lwd=3, col="indianred4")
mean(a)
length(which(a==0))
length(which(a<20))

