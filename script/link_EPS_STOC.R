rm(list=ls())

# Packages necessaires
devtools::install_deps(upgrade="never")

# Load fonctions importantes
devtools::load_all() 

#EPS
load(here::here("output","CLC_EPS.RData"))

# STOC
load(here::here("output","CLC_STOC.RData"))


# Creation des couche de points
dsf_STOC <-  Mesanges::give_point(CLC_STOC)
dsf_EPS <-  Mesanges::give_point(CLC_EPS)


# Point par par station STOC
EPS_by_STOC <- list()
for (i in 1:nrow(CLC_STOC)) {
  EPS_by_STOC[[i]] <- Mesanges::give_point_by_site(dsf_STOC[i,], 25, CLC_EPS, dsf_EPS)
}

save(EPS_by_STOC, file  = here::here("output","EPS_by_STOC.RData"))


###### Carte

# Load Corine Land Cover
shp_CLC <- sf::st_read(dsn   = here::here("data","CLC12_FR_RGF_SHP"),
                       layer = "CLC12_FR_RGF") %>%
  sf::st_transform(shp,crs = 2154)

load(here::here("output","EPS_by_STOC.RData"))

Mesanges::plot_carte(25, dsf_STOC, CLC_EPS, 25, shp_CLC)

##### Info en supplément

EPS_by_STOC <- list()
for (i in 1:nrow(CLC_STOC)) {
  EPS_by_STOC[[i]] <- Mesanges::give_point_by_site(dsf_STOC[i,], 25, CLC_EPS)
}

a <- NULL
for (i in 1:nrow(CLC_STOC)) {
  a[i] <- length(EPS_by_STOC[[i]])
}

hist(a, breaks= 50, col= "ivory3", xlab="Nombre de point d'écoute", ylab="Fréquence", main = "Nombre de points d'écoute dans un buffer de 20 km autour de la station STOC")
abline(v=20, lty=2, lwd=3, col="indianred4")
mean(a)
length(which(a==0))
length(which(a<20))

