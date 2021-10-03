#creation du RC
rrtools::use_compendium("../Fauvettes/", open = FALSE)

#Creation de mon fichier de fonction
usethis::use_r(name="function.R") #cree dossier R avec le fichier fonction R
#Creation des dossiers
dir.create("data")
dir.create("script")
dir.create("output")
usethis::use_git_ignore("data/")

#Ajout de packages utiles
usethis::use_package("devtools")
usethis::use_package("sp")
usethis::use_package("ggplot2")
usethis::use_package("here")
usethis::use_package("AER")
usethis::use_package("nimble")
usethis::use_pipe()
devtools::document()

#verifier notre package
devtools::check()

# readme ?
# rrtools::use_readme_rmd()
