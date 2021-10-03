##########################
#   Index fauvettes      #
##########################
rm(list=ls())
color<-c("#FF8830","#A6B06D","#589482","#8C2423")

# Packages necessaires
devtools::install_deps(upgrade="never")

# Load fonctions importantes
devtools::load_all() 

# Les donnees EPS
load(here::here("output","data_EPS.RData"))
data_EPS <- data_EPS %>%
  dplyr::filter(data_EPS$annee != 1914) %>%
  dplyr::rename(
    long = 'longitude_wgs84',
    lat  = 'latitude_wgs84') %>%
  dplyr::filter(!is.na(long)) %>%
  dplyr::filter(!is.na(lat)) %>%
  replace(is.na(.), 0)

# Type d'habitat pour chaque station STOC
load(here::here("output","CLC_STOC.RData"))

# Type d'habitat pour chaque station EPS
load(here::here("output","CLC_EPS.RData"))

# Les points d'ecoute (juste leurs identifiants) pour chaque station STOC dans un rayon de 25km
load(here::here("output","EPS_by_STOC.RData"))

# On lie les tables CLC_EPS et data_EPS
EPS <- dplyr::left_join(data_EPS, CLC_EPS, by = c("point", "lat", "long")) %>%
  dplyr::rename(CLC_EPS = 'CLC_hab')  %>%
  dplyr::select(-CLC_1) %>%
  dplyr::select(-CLC_2)

# Table finale
# selectionner les point concernant la station,
# ajouter une colonne avec le nom de la station
# lier les nouveaux tableaux

data <- data.frame("carre" = NA, "annee" = NA, "point" = NA, long= as.numeric("NA"), lat =as.numeric("NA"), PARCAE = as.numeric("NA"), PARMAJ = as.numeric("NA"), SYLATR = as.numeric("NA"), SYLBOR = as.numeric("NA"), CLC_EPS = as.numeric("NA"), ID_PROG = "NA", CLC_STOC = "NA")
# warnings "NA" pas importants
for (i in 1:nrow(CLC_STOC)) {
  data_add <- EPS %>%
    dplyr::filter(EPS$point %in% EPS_by_STOC[[i]]) %>%
    tibble::add_column(ID_PROG = as.character(CLC_STOC$ID_PROG[i])) %>%
    tibble::add_column(CLC_STOC = as.character(CLC_STOC$CLC_hab[i]))
  
  data <- dplyr::bind_rows(data, data_add)
}

rm(data_add, EPS_by_STOC, CLC_EPS, CLC_STOC, EPS, i)

# Avec les donnees de TITS
#############################

load(here::here("output","hvie_ID_PROG_syl.RData"))
ID_PROG_syl <- hvie_ID_PROG_syl$ID_PROG
# On ne garde que les sites utilises dans les donnees parmaj
data_syl <- data %>%
  dplyr::filter(data$ID_PROG %in% ID_PROG_syl)

rm(hvie_ID_PROG_syl, ID_PROG_syl, data_EPS)

# On supprime les points qui n'ont jamais contacté de mesanges
data_syl <- Mesanges::supr_point(data_syl, "SYL")
# on supprime les doublons
data_syl <- data_syl %>% dplyr::distinct(carre,annee,point,.keep_all = TRUE)

save(data_syl,file=here::here('output','data_syl.RData'))

##############################################################################################
# CREATION INDEX
##############################################################################################

### SYLBOR 
index_sylbor_point <- Mesanges::index_new(data_syl, "SYLBOR")
save(index_sylbor_point,file=here::here("output","index_sylbor_point.RData"))
load(here::here("output","index_sylbor_point.RData"))
Mesanges::plot_index(index_sylbor_point[[1]], sqrt(index_sylbor_point[[2]]), color[3], "sylbor", "point" )

### SYLATR
index_sylatr_point <- Mesanges::index_new(data_syl, "SYLATR")
save(index_sylatr_point,file=here::here("output","index_sylatr_point.RData"))
load(here::here("output","index_sylatr_point.RData"))
Mesanges::plot_index(index_sylatr_point[[1]], sqrt(index_sylatr_point[[2]]), color[1], "sylatr", "point" )


##############################################################################################
# PLOT INDEX
##############################################################################################

load(here::here("output","index_sylbor_point.RData"))
ind_gw <- list()
ind_gw$gammat_gw <- t(apply(index_sylbor_point[[3]],c(2,3), mean, na.rm = TRUE))
ind_gw$l_025_gw  <- apply(ind_gw$gammat_gw,2,quantile,probs=0.025)
ind_gw$l_25_gw   <- apply(ind_gw$gammat_gw,2,quantile,probs=0.25)
ind_gw$l_50_gw   <- apply(ind_gw$gammat_gw,2,quantile,probs=0.50)
ind_gw$l_75_gw   <- apply(ind_gw$gammat_gw,2,quantile,probs=0.75)
ind_gw$l_975_gw  <- apply(ind_gw$gammat_gw,2,quantile,probs=0.975)

mean(ind_gw$gammat_gw)
sd(ind_gw$gammat_gw)

load(here::here("output","index_sylatr_point.RData"))
ind_bc <- list()
ind_bc$gammat_bc <- t(apply(index_sylatr_point[[3]],c(2,3), mean, na.rm = TRUE))
ind_bc$l_025_bc  <- apply(ind_bc$gammat_bc,2,quantile,probs=0.025)
ind_bc$l_25_bc   <- apply(ind_bc$gammat_bc,2,quantile,probs=0.25)
ind_bc$l_50_bc   <- apply(ind_bc$gammat_bc,2,quantile,probs=0.50)
ind_bc$l_75_bc   <- apply(ind_bc$gammat_bc,2,quantile,probs=0.75)
ind_bc$l_975_bc  <- apply(ind_bc$gammat_bc,2,quantile,probs=0.975)

mean(ind_bc$gammat_bc)
sd(ind_bc$gammat_bc)

color<-c("#2B2D42","#E4BB97","#D08C60","indianred4")
color_gw <- color[3]
color_bc <- color[1]

color.transparent_gw <- adjustcolor(color_gw, alpha.f = 0.4)
color.transparent_2_gw<- adjustcolor(color_gw, alpha.f = 0.2)
color.transparent_bc <- adjustcolor(color_bc, alpha.f = 0.4)
color.transparent_2_bc<- adjustcolor(color_bc, alpha.f = 0.2)


K=18
par(mfrow=c(1,1))
par(cex=0.9, mai=c(0.7,0.8,0.5,0.1))

plot(ind_gw$l_50_gw,type='l',axes=F,lwd =2,col="ivory4",ylim=c(0,2.5),ylab="Mean number of birds per count point",xlab="",main="")
xx<- c(1:K,K:1)
yy <- c(ind_gw$l_025_gw[1:K],ind_gw$l_975_gw[K:1])
polygon(xx,yy,col=color.transparent_2_gw, border=NA)
xx<- c(1:K,K:1)
yy <- c(ind_bc$l_025_bc[1:K],ind_bc$l_975_bc[K:1])
polygon(xx,yy,col=color.transparent_2_bc, border=NA)
xx<- c(1:K,K:1)
yy <- c(ind_gw$l_25_gw[1:K],ind_gw$l_75_gw[K:1])
polygon(xx,yy,col=color.transparent_gw, border=NA)
xx<- c(1:K,K:1)
yy <- c(ind_bc$l_25_bc[1:K],ind_bc$l_75_bc[K:1])
polygon(xx,yy,col=color.transparent_bc, border=NA)
lines(ind_gw$l_50_gw,lwd =2, col=color_gw)
lines(ind_bc$l_50_bc,lwd =2, col=color_bc)
axis(1, at=seq(1,K,4),labels=c(seq(2002,2019,4)))
axis(side =2, cex.axis=1, las=2)
legend(13,2.67,legend=c("Blackcap", "Garden warbler"), lwd=2,lty=c(1,1),col=c(color[c(1,3)]), bty='n')


