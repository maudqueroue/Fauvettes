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

# type d'habitat
#--------------------
# On ajoute les infos contenu dans CLC_STOC aux hvies 
load(here::here("output","CLC_STOC.RData"))
hvie_ID_PROG_syl <- dplyr::left_join(hvie_ID_PROG_syl, CLC_STOC, by = c("ID_PROG")) 

rm(CLC_STOC)  


# 11. Creation vecteur pour lier hvie individus et hvie site
#---------
hvie_sylbor <- Mesanges::link_hvie(hvie_sylbor, hvie_ID_PROG_syl)
hvie_sylatr <- Mesanges::link_hvie(hvie_sylatr, hvie_ID_PROG_syl)

save(hvie_ID_PROG_syl,file=here::here('output',"hvie_ID_PROG_syl.RData"))
save(hvie_sylbor,file=here::here('output',"hvie_sylbor.RData"))
save(hvie_sylatr,file=here::here('output',"hvie_sylatr.RData"))

