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

#STOC
#-----------------
# Load sites
load(here::here('output','hvie_ID_PROG_syl.RData'))

data_STOC <-  hvie_ID_PROG_syl %>%
  dplyr::select('lat','long')
data_STOC$type <- "STOC"

# EPS
#--------------
load(here::here('output','data_syl.RData'))

data_syl <- data_syl %>%
  dplyr::distinct(point, .keep_all=TRUE)

data_EPS <- data_syl %>%
  dplyr::select('lat','long')

data_EPS$type <- "EPS"

#data complet 
#---------------
data <- rbind(data_EPS,data_STOC)

# Carte
#----------------
ggplot2::ggplot(Fr) +
  ggplot2::xlab("") + 
  ggplot2::ylab("") +
  ggplot2::xlim(-4.5, 7.8) +
  ggplot2::ylim(42.5,51) +
  ggplot2::geom_sf() +
  ggplot2::geom_segment(ggplot2::aes(x = -2, y = 45, xend = 7.5, yend = 45),linetype="dashed",size=1.3) +
  ggplot2::geom_point(data = data, ggplot2::aes(x= long, y=lat, color=type, shape = type, size=type)) +
  ggplot2::scale_color_manual(values=c('#C75146','#540804'),name='',labels=c('Count points','CR sites')) +
  ggplot2::scale_shape_manual(values=c(19,17,17),name='',labels=c('Count points','CR sites')) +
  ggplot2::scale_size_manual(values=c(0.9,3,3),name='',labels=c('Count points','CR sites')) +
  ggplot2::geom_text(x=-2.8, y=45, label="45°N",size=3)+
  ggplot2::theme_minimal()


# Informations donnees tits
#-----------------------------------
########## EPS

load(here::here('output','data_syl.RData'))
length(unique(data_syl$point))
sum(data_syl$SYLATR)
sum(data_syl$SYLBOR)

EPS_gw <- NULL
EPS_bc <- NULL

year <- seq(2001,2019,1)
for(i in 1 :19){
  sub_point_EPS <- subset(data_syl,data_syl$annee==year[i])
  EPS_gw[i] <- sum(sub_point_EPS$SYLBOR)/nrow(sub_point_EPS)
  EPS_bc[i] <- sum(sub_point_EPS$SYLATR)/nrow(sub_point_EPS)
  rm(sub_point_EPS)
}

mean(EPS_gw)
mean(EPS_bc)

# STOC
load(here::here('output','hvie_sylbor.RData'))
load(here::here('output','hvie_sylatr.RData'))

#Nb individus
length(unique(hvie_sylbor$ID))
length(unique(hvie_sylatr$ID))


# GW
#-------------------
#Infos
temp <- 1:19
AGE_gw <- NULL
nbrecapt_gw <- NULL

for (i in 1:nrow(hvie_sylbor)){
  AGE_gw[i] <- hvie_sylbor[i,min(temp[hvie_sylbor[i,1:19]>0])]
  nbrecapt_gw[i] <- length(which(hvie_sylbor[i,1:19]>0))-1
}

# Calcul nb recapture
length(which(nbrecapt_gw>0))
# info complet
info_gw <- cbind(AGE_gw,nbrecapt_gw)
# Calcul nb jeunes recaptured
length(which(info_gw[,1]=="1" & info_gw[,2]>0))
# Calcul nb adults recaptured
sum(nbrecapt_gw) - length(which(info_gw[,1]=="1" & info_gw[,2]>0))

# BC
#-------------------
#Infos
temp <- 1:19
AGE_bc <- NULL
nbrecapt_bc <- NULL

for (i in 1:nrow(hvie_sylatr)){
  AGE_bc[i] <- hvie_sylatr[i,min(temp[hvie_sylatr[i,1:19]>0])]
  nbrecapt_bc[i] <- length(which(hvie_sylatr[i,1:19]>0))-1
}

# Calcul nb recapture
length(which(nbrecapt_bc>0))
# info complet
info_bc <- cbind(AGE_bc,nbrecapt_bc)
# Calcul nb jeunes recaptured
length(which(info_bc[,1]=="1" & info_bc[,2]>0))
# Calcul nb adults recaptured
sum(nbrecapt_bc) - length(which(info_bc[,1]=="1" & info_bc[,2]>0))

# Mean nb capt STOC
#--------------------
load(here::here('output','hvie_ID_PROG_syl.RData'))
hvie_site <- as.matrix(hvie_ID_PROG_syl[1:19])
hvie_site[which(hvie_site>0)] <- 1

nb_capt_site <- apply(hvie_site[,1:19],2,sum)

sylbor <- as.matrix(hvie_sylbor[,1:19])
sylbor <- apply(sylbor,2,as.numeric)
sylbor[which(sylbor>0)] <- 1
nb_capt_sylbor <- apply(sylbor,2,sum)
mean(nb_capt_sylbor/nb_capt_site)

sylatr <- as.matrix(hvie_sylatr[,1:19])
sylatr <- apply(sylatr,2,as.numeric)
sylatr[which(sylatr>0)] <- 1
nb_capt_sylatr <- apply(sylatr,2,sum)
mean(nb_capt_sylatr/nb_capt_site)

