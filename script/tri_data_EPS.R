rm(list=ls())

# Packages necessaires
devtools::install_deps(upgrade="never")

# Load fonctions importantes
devtools::load_all() 

# Les data
data <- read.table(here::here("data","STOC_EPS.csv"),head=T,sep=";") 

data_geo <- data %>%
        dplyr::select("point","carre","longitude_wgs84","latitude_wgs84") %>%
          dplyr::distinct(point, .keep_all = T)

data_sp <- data %>%
             dplyr::select("point","carre","annee","code_sp","abondance") %>%
             tidyr::pivot_wider(names_from = code_sp, values_from = abondance) 
             
# Les data sur les carres
carre <- read.table(here::here("data","EPS_carres.csv"),head=T,sep=";")

carre <- carre %>%
  dplyr::filter(carre$annee != 1914)
unique(carre$annee)

id_carre_tot <- carre %>%
                      dplyr::distinct(id_carre)

# On reformate les data dans new pour ajouter les zero 
new <- data.frame("carre" = NA, "annee" = NA, "point" = NA)


for(i in 1:nrow(id_carre_tot)) {
  
  # nom carre 
  name_carre <- id_carre_tot$id_carre[i]
  
  # Nombre de point dans le carre
  data_carre <- data %>%
    dplyr::filter(carre %in% name_carre)
  name_point <- dplyr::distinct(as.data.frame(data_carre$point)) %>%
    dplyr::rename(name_point = 'data_carre$point')
  nb_point <- dplyr::n_distinct(name_point)
  
  # Nombre d annees effectuees
  annee_carre <- carre %>%
    dplyr::filter(id_carre %in% name_carre)
  name_annee <- dplyr::distinct(as.data.frame(annee_carre$annee)) %>%
    dplyr::rename(name_annee = 'annee_carre$annee')
  nb_annee <- dplyr::n_distinct(name_annee)
 
  if(nb_point > 0) {
  x <- data.frame("carre" = rep(name_carre,nb_point*nb_annee), "annee" = rep(name_annee$name_annee,each = nb_point), "point" = rep(name_point$name_point,nb_annee))
   }
  
  if(nb_point == 0) {
    x <- data.frame("carre" = rep(name_carre,nb_annee), "annee" = name_annee$name_annee, "point" = rep(NA,nb_annee))
  }
  
  new <- dplyr::bind_rows(new,x)
  rm(x, nb_annee, name_annee, nb_point, name_point, name_carre, data_carre, annee_carre)
} 

# On fait le lien entre les point et leur geolocalisation
new_geo <- dplyr::left_join(new, data_geo, by = c("carre", "point"))
# On ajoute les comptage de sp
data_EPS <- dplyr::left_join(new_geo, data_sp, by =c ("carre", "point", "annee")) 

rm(data_geo, data, carre, data_sp, id_carre_tot, new, new_geo, i)

# On sauve le fichier
save(data_EPS, file = here::here("output","data_EPS.RData"))
